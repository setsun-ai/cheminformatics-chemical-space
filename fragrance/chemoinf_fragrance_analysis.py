#!/usr/bin/env python3
"""
Chemoinformatics project pipeline for fragrance/flavour molecules.

Input CSV columns required:
    id, group, compound_name, smiles
Optional columns:
    odor_note, canonical_smiles_rdkit, etc.

Outputs:
    - validated molecule table
    - RDKit physicochemical descriptors
    - MACCS and Morgan fingerprints
    - PCA / t-SNE / optional UMAP coordinates and plots
    - hierarchical clustering dendrograms and cluster labels
    - simple quantitative comparison tables
    - 2D molecule grids
    - 3D SDF file for Jmol / Avogadro / PyMOL-style viewers

Run example:
    python chemoinf_fragrance_analysis.py \
        --input chemoinf_fragrance_dataset_v1.csv \
        --outdir chemoinf_results

Written with ChatGPT assistance.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist, squareform

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.preprocessing import StandardScaler

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Draw, Lipinski, MACCSkeys, rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

try:
    import umap  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    umap = None

warnings.filterwarnings("ignore", category=UserWarning)


PHYS_DESC_COLUMNS = [
    "MolWt",
    "ExactMolWt",
    "MolLogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotatableBonds",
    "RingCount",
    "AromaticRings",
    "AliphaticRings",
    "HeavyAtomCount",
    "FractionCSP3",
    "NumHeteroatoms",
    "NumValenceElectrons",
    "NHOHCount",
    "NOCount",
]


def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for sub in [
        "tables",
        "plots_pca",
        "plots_tsne",
        "plots_umap",
        "plots_hclust",
        "structures_2d",
        "structures_3d",
    ]:
        (path / sub).mkdir(exist_ok=True)


def load_dataset(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path)
    required = {"id", "group", "compound_name", "smiles"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    df = df.copy()
    df["id"] = df["id"].astype(str)
    df["group"] = df["group"].astype(str)
    df["compound_name"] = df["compound_name"].astype(str)
    df["smiles"] = df["smiles"].astype(str)
    return df


def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def validate_molecules(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[Chem.Mol]]:
    mols: List[Optional[Chem.Mol]] = []
    canonical = []
    valid = []

    for smi in df["smiles"]:
        mol = mol_from_smiles(smi)
        mols.append(mol)
        valid.append(mol is not None)
        canonical.append(Chem.MolToSmiles(mol, canonical=True) if mol is not None else np.nan)

    out = df.copy()
    out["rdkit_valid"] = valid
    out["canonical_smiles"] = canonical

    valid_df = out[out["rdkit_valid"]].reset_index(drop=True)
    valid_mols = [m for m in mols if m is not None]

    duplicate_count = valid_df["canonical_smiles"].duplicated().sum()
    if duplicate_count > 0:
        print(f"WARNING: {duplicate_count} duplicated canonical SMILES found.", file=sys.stderr)

    return valid_df, valid_mols


def calc_physchem_descriptors(mols: List[Chem.Mol]) -> pd.DataFrame:
    rows = []
    for mol in mols:
        rows.append(
            {
                "MolWt": Descriptors.MolWt(mol),
                "ExactMolWt": Descriptors.ExactMolWt(mol),
                "MolLogP": Descriptors.MolLogP(mol),
                "TPSA": rdMolDescriptors.CalcTPSA(mol),
                "HBD": Lipinski.NumHDonors(mol),
                "HBA": Lipinski.NumHAcceptors(mol),
                "RotatableBonds": Lipinski.NumRotatableBonds(mol),
                "RingCount": Lipinski.RingCount(mol),
                "AromaticRings": Lipinski.NumAromaticRings(mol),
                "AliphaticRings": Lipinski.NumAliphaticRings(mol),
                "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
                "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
                "NumHeteroatoms": Lipinski.NumHeteroatoms(mol),
                "NumValenceElectrons": Descriptors.NumValenceElectrons(mol),
                "NHOHCount": Lipinski.NHOHCount(mol),
                "NOCount": Lipinski.NOCount(mol),
            }
        )
    return pd.DataFrame(rows)


def bitvect_to_numpy(fp, n_bits: int) -> np.ndarray:
    arr = np.zeros((n_bits,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def calc_maccs_fingerprints(mols: List[Chem.Mol]) -> pd.DataFrame:
    rows = []
    for mol in mols:
        fp = MACCSkeys.GenMACCSKeys(mol)  # usually 167 bits, bit 0 unused
        arr = bitvect_to_numpy(fp, fp.GetNumBits())
        rows.append(arr)
    cols = [f"MACCS_{i:03d}" for i in range(len(rows[0]))]
    return pd.DataFrame(rows, columns=cols)


def calc_morgan_fingerprints(mols: List[Chem.Mol], radius: int = 2, n_bits: int = 2048) -> pd.DataFrame:
    gen = GetMorganGenerator(radius=radius, fpSize=n_bits)
    rows = []
    for mol in mols:
        fp = gen.GetFingerprint(mol)
        rows.append(bitvect_to_numpy(fp, n_bits))
    cols = [f"Morgan_r{radius}_{i:04d}" for i in range(n_bits)]
    return pd.DataFrame(rows, columns=cols)


def save_table_with_id(df_meta: pd.DataFrame, features: pd.DataFrame, out_path: Path) -> None:
    keep = [c for c in ["id", "group", "compound_name", "odor_note", "smiles", "canonical_smiles"] if c in df_meta.columns]
    table = pd.concat([df_meta[keep].reset_index(drop=True), features.reset_index(drop=True)], axis=1)
    table.to_csv(out_path, index=False)


def plot_embedding(
    coords: np.ndarray,
    df_meta: pd.DataFrame,
    title: str,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    explained: Optional[Tuple[float, float]] = None,
) -> None:
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    groups = list(pd.Series(df_meta["group"]).dropna().unique())
    for group in groups:
        mask = (df_meta["group"].values == group)
        ax.scatter(coords[mask, 0], coords[mask, 1], label=group, alpha=0.82, s=45)

    if explained is not None:
        xlabel = f"{xlabel} ({explained[0] * 100:.1f}% var.)"
        ylabel = f"{ylabel} ({explained[1] * 100:.1f}% var.)"

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=True)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def run_pca(
    X: np.ndarray,
    df_meta: pd.DataFrame,
    name: str,
    outdir: Path,
    scale: bool = False,
) -> Tuple[pd.DataFrame, Tuple[float, float]]:
    X_work = X.astype(float)
    if scale:
        X_work = StandardScaler().fit_transform(X_work)

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_work)
    evr = tuple(pca.explained_variance_ratio_[:2])

    coord_df = df_meta[["id", "group", "compound_name"]].copy()
    coord_df[f"{name}_PCA1"] = coords[:, 0]
    coord_df[f"{name}_PCA2"] = coords[:, 1]
    coord_df.to_csv(outdir / "tables" / f"coords_pca_{name}.csv", index=False)

    plot_embedding(
        coords,
        df_meta,
        title=f"PCA — {name}",
        xlabel="PC1",
        ylabel="PC2",
        out_path=outdir / "plots_pca" / f"pca_{name}.png",
        explained=evr,
    )
    return coord_df, evr


def jaccard_distance_matrix(X_bits: np.ndarray) -> np.ndarray:
    return squareform(pdist(X_bits.astype(bool), metric="jaccard"))


def run_tsne(
    X: np.ndarray,
    df_meta: pd.DataFrame,
    name: str,
    outdir: Path,
    metric: str,
    perplexities: Iterable[int] = (5, 15, 30),
    scale: bool = False,
    max_iter: int = 750,
) -> List[pd.DataFrame]:
    outputs = []
    n = X.shape[0]

    if metric == "precomputed_jaccard":
        D = jaccard_distance_matrix(X)
        X_input = D
        metric_arg = "precomputed"
        init_arg = "random"
    else:
        X_input = X.astype(float)
        if scale:
            X_input = StandardScaler().fit_transform(X_input)
        metric_arg = metric
        init_arg = "pca"

    for perp in perplexities:
        if perp >= n:
            continue
        tsne = TSNE(
            n_components=2,
            perplexity=perp,
            metric=metric_arg,
            init=init_arg,
            learning_rate="auto",
            max_iter=max_iter,
            random_state=42,
        )
        coords = tsne.fit_transform(X_input)

        coord_df = df_meta[["id", "group", "compound_name"]].copy()
        coord_df[f"{name}_TSNE1_p{perp}"] = coords[:, 0]
        coord_df[f"{name}_TSNE2_p{perp}"] = coords[:, 1]
        coord_df.to_csv(outdir / "tables" / f"coords_tsne_{name}_perplexity_{perp}.csv", index=False)

        plot_embedding(
            coords,
            df_meta,
            title=f"t-SNE — {name}, perplexity={perp}",
            xlabel="t-SNE 1",
            ylabel="t-SNE 2",
            out_path=outdir / "plots_tsne" / f"tsne_{name}_p{perp}.png",
        )
        outputs.append(coord_df)
    return outputs


def run_umap_optional(
    X: np.ndarray,
    df_meta: pd.DataFrame,
    name: str,
    outdir: Path,
    metric: str,
    n_neighbors_values: Iterable[int] = (5, 15, 30),
    min_dist_values: Iterable[float] = (0.1, 0.5),
    scale: bool = False,
) -> List[pd.DataFrame]:
    if umap is None:
        return []

    X_input = X.astype(float)
    if scale:
        X_input = StandardScaler().fit_transform(X_input)

    outputs = []
    n = X.shape[0]
    for nn in n_neighbors_values:
        if nn >= n:
            continue
        for min_dist in min_dist_values:
            reducer = umap.UMAP(
                n_components=2,
                n_neighbors=nn,
                min_dist=min_dist,
                metric=metric,
                random_state=42,
            )
            coords = reducer.fit_transform(X_input)
            suffix = f"nn{nn}_mindist{str(min_dist).replace('.', '_')}"

            coord_df = df_meta[["id", "group", "compound_name"]].copy()
            coord_df[f"{name}_UMAP1_{suffix}"] = coords[:, 0]
            coord_df[f"{name}_UMAP2_{suffix}"] = coords[:, 1]
            coord_df.to_csv(outdir / "tables" / f"coords_umap_{name}_{suffix}.csv", index=False)

            plot_embedding(
                coords,
                df_meta,
                title=f"UMAP — {name}, n_neighbors={nn}, min_dist={min_dist}",
                xlabel="UMAP 1",
                ylabel="UMAP 2",
                out_path=outdir / "plots_umap" / f"umap_{name}_{suffix}.png",
            )
            outputs.append(coord_df)
    return outputs


def run_hierarchical_clustering(
    X: np.ndarray,
    df_meta: pd.DataFrame,
    name: str,
    outdir: Path,
    kind: str,
    scale: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    labels_true = pd.Categorical(df_meta["group"]).codes

    if kind == "physchem":
        X_work = X.astype(float)
        if scale:
            X_work = StandardScaler().fit_transform(X_work)
        Z = linkage(X_work, method="ward", metric="euclidean")
        cluster_labels = fcluster(Z, t=3, criterion="maxclust")
        sil_true = silhouette_score(X_work, labels_true, metric="euclidean")
        sil_hclust = silhouette_score(X_work, cluster_labels, metric="euclidean")
    elif kind == "fingerprint":
        condensed = pdist(X.astype(bool), metric="jaccard")
        D = squareform(condensed)
        Z = linkage(condensed, method="average")
        cluster_labels = fcluster(Z, t=3, criterion="maxclust")
        sil_true = silhouette_score(D, labels_true, metric="precomputed")
        sil_hclust = silhouette_score(D, cluster_labels, metric="precomputed")
    else:
        raise ValueError("kind must be either 'physchem' or 'fingerprint'")

    ari = adjusted_rand_score(labels_true, cluster_labels)
    nmi = normalized_mutual_info_score(labels_true, cluster_labels)

    cluster_df = df_meta[["id", "group", "compound_name"]].copy()
    cluster_df[f"hclust_{name}_cluster_k3"] = cluster_labels
    cluster_df.to_csv(outdir / "tables" / f"hclust_clusters_{name}.csv", index=False)

    labels = [f"{r.id} | {r.group} | {r.compound_name}" for r in df_meta.itertuples(index=False)]
    fig = plt.figure(figsize=(13, 7))
    ax = fig.add_subplot(111)
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=6, ax=ax)
    ax.set_title(f"Hierarchical clustering — {name}")
    ax.set_ylabel("Distance")
    fig.tight_layout()
    fig.savefig(outdir / "plots_hclust" / f"hclust_dendrogram_{name}.png", dpi=300)
    plt.close(fig)

    metrics = {
        "ARI_HCA_k3_vs_group": float(ari),
        "NMI_HCA_k3_vs_group": float(nmi),
        "silhouette_true_groups": float(sil_true),
        "silhouette_HCA_k3": float(sil_hclust),
    }
    return cluster_df, metrics


def make_group_descriptor_means(df_meta: pd.DataFrame, physchem: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    table = pd.concat([df_meta[["id", "group", "compound_name"]], physchem], axis=1)
    means = table.groupby("group")[PHYS_DESC_COLUMNS].mean().round(3)
    medians = table.groupby("group")[PHYS_DESC_COLUMNS].median().round(3)
    means.to_csv(outdir / "tables" / "physchem_group_means.csv")
    medians.to_csv(outdir / "tables" / "physchem_group_medians.csv")
    return means


def draw_molecule_grids(df_meta: pd.DataFrame, mols: List[Chem.Mol], outdir: Path) -> None:
    grid_dir = outdir / "structures_2d"
    legends_all = [f"{r.id}\n{r.compound_name}" for r in df_meta.itertuples(index=False)]

    img = Draw.MolsToGridImage(mols, legends=legends_all, molsPerRow=5, subImgSize=(260, 190), useSVG=False)
    img.save(grid_dir / "all_molecules_grid.png")

    for group in df_meta["group"].unique():
        mask = df_meta["group"].values == group
        group_mols = [mol for mol, keep in zip(mols, mask) if keep]
        group_meta = df_meta.loc[mask].reset_index(drop=True)
        legends = [f"{r.id}\n{r.compound_name}" for r in group_meta.itertuples(index=False)]
        img = Draw.MolsToGridImage(group_mols, legends=legends, molsPerRow=5, subImgSize=(260, 190), useSVG=False)
        img.save(grid_dir / f"{group}_molecules_grid.png")

    # Smaller representative grid for slides: first 6 molecules from each group.
    rep_mols = []
    rep_legends = []
    for group in df_meta["group"].unique():
        idxs = list(np.where(df_meta["group"].values == group)[0])[:6]
        for i in idxs:
            rep_mols.append(mols[i])
            rep_legends.append(f"{df_meta.loc[i, 'group']}\n{df_meta.loc[i, 'compound_name']}")
    img = Draw.MolsToGridImage(rep_mols, legends=rep_legends, molsPerRow=6, subImgSize=(240, 180), useSVG=False)
    img.save(grid_dir / "representative_molecules_grid.png")


def make_3d_conformer(mol: Chem.Mol, random_seed: int = 42) -> Optional[Chem.Mol]:
    mol3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    if hasattr(params, "maxAttempts"):
        params.maxAttempts = 1000
    elif hasattr(params, "maxIterations"):
        params.maxIterations = 1000
    try:
        status = AllChem.EmbedMolecule(mol3d, params)
        if status != 0:
            return None
        try:
            AllChem.UFFOptimizeMolecule(mol3d, maxIters=500)
        except Exception:
            # Geometry is still usable for visual inspection even if UFF optimization fails.
            pass
        return mol3d
    except Exception:
        return None


def write_3d_sdf(df_meta: pd.DataFrame, mols: List[Chem.Mol], outdir: Path) -> pd.DataFrame:
    sdf_path = outdir / "structures_3d" / "fragrance_molecules_3d_for_jmol.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    rows = []
    for i, (row, mol) in enumerate(zip(df_meta.itertuples(index=False), mols)):
        mol3d = make_3d_conformer(mol, random_seed=42 + i)
        success = mol3d is not None
        if mol3d is not None:
            mol3d.SetProp("_Name", f"{row.id}_{row.compound_name}")
            mol3d.SetProp("id", str(row.id))
            mol3d.SetProp("group", str(row.group))
            mol3d.SetProp("compound_name", str(row.compound_name))
            writer.write(mol3d)
        rows.append({"id": row.id, "group": row.group, "compound_name": row.compound_name, "3d_conformer_success": success})
    writer.close()

    status_df = pd.DataFrame(rows)
    status_df.to_csv(outdir / "tables" / "sdf_3d_generation_status.csv", index=False)
    return status_df


def write_summary_report(
    df_meta: pd.DataFrame,
    outdir: Path,
    evr_table: pd.DataFrame,
    metrics_table: pd.DataFrame,
    umap_available: bool,
    sdf_status: Optional[pd.DataFrame],
) -> None:
    counts = df_meta["group"].value_counts().sort_index()
    duplicates = int(df_meta["canonical_smiles"].duplicated().sum()) if "canonical_smiles" in df_meta else 0

    lines = []
    lines.append("# Chemoinformatics fragrance analysis — automatic report\n")
    lines.append("## Dataset\n")
    lines.append(f"- Valid molecules: **{len(df_meta)}**")
    lines.append(f"- Duplicated canonical SMILES: **{duplicates}**")
    lines.append("- Group counts:")
    for group, count in counts.items():
        lines.append(f"  - {group}: {count}")

    lines.append("\n## Descriptor sets\n")
    lines.append("- Physchem: scaled RDKit physicochemical descriptors for PCA/t-SNE/HCA.")
    lines.append("- MACCS: 167-bit structural keys, analysed mainly with Jaccard distance.")
    lines.append("- Morgan: radius 2, 2048-bit circular fingerprints, analysed mainly with Jaccard distance.")

    lines.append("\n## PCA explained variance\n")
    lines.append(evr_table.to_markdown(index=False))

    lines.append("\n## Clustering / separability metrics\n")
    lines.append(metrics_table.to_markdown(index=False))

    lines.append("\n## UMAP\n")
    if umap_available:
        lines.append("UMAP plots were generated because `umap-learn` is installed.")
    else:
        lines.append("UMAP was skipped because `umap-learn` is not installed. Install it with `pip install umap-learn` or `conda install -c conda-forge umap-learn`.")

    lines.append("\n## 3D SDF for Jmol\n")
    if sdf_status is not None:
        success_count = int(sdf_status["3d_conformer_success"].sum())
        lines.append(f"Generated 3D conformers for **{success_count}/{len(sdf_status)}** molecules.")
        lines.append("Open `structures_3d/fragrance_molecules_3d_for_jmol.sdf` in Jmol/Avogadro for 3D inspection.")
    else:
        lines.append("3D SDF generation was skipped.")

    lines.append("\n## Recommended first plots for presentation\n")
    lines.append("1. `structures_2d/representative_molecules_grid.png`")
    lines.append("2. `plots_pca/pca_physchem.png`")
    lines.append("3. `plots_pca/pca_maccs.png`")
    lines.append("4. `plots_pca/pca_morgan_r2_2048.png`")
    lines.append("5. one t-SNE plot for Morgan, preferably perplexity 15 or 30")
    lines.append("6. `plots_hclust/hclust_dendrogram_morgan_r2_2048.png`")

    (outdir / "analysis_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_int_list(text: str) -> List[int]:
    values = []
    for item in text.split(","):
        item = item.strip()
        if item:
            values.append(int(item))
    if not values:
        raise ValueError("At least one integer value is required.")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description="Chemoinformatics fragrance/flavour molecule analysis pipeline.")
    parser.add_argument("--input", required=True, help="Input CSV with id, group, compound_name, smiles columns.")
    parser.add_argument("--outdir", default="chemoinf_results", help="Output directory.")
    parser.add_argument("--morgan-bits", type=int, default=2048, help="Morgan fingerprint size.")
    parser.add_argument("--morgan-radius", type=int, default=2, help="Morgan fingerprint radius.")
    parser.add_argument("--skip-images", action="store_true", help="Skip 2D molecule grid generation.")
    parser.add_argument("--skip-3d", action="store_true", help="Skip 3D SDF generation for Jmol.")
    parser.add_argument("--skip-tsne", action="store_true", help="Skip t-SNE plots. Useful for a quick first run.")
    parser.add_argument("--tsne-perplexities", default="5,15,30", help="Comma-separated t-SNE perplexities, e.g. 5,15,30.")
    parser.add_argument("--tsne-iterations", type=int, default=750, help="Number of t-SNE optimization iterations.")
    args = parser.parse_args()
    tsne_perplexities = parse_int_list(args.tsne_perplexities)

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    print(f"Loading dataset: {input_path}")
    raw_df = load_dataset(input_path)
    df_meta, mols = validate_molecules(raw_df)

    invalid_df = raw_df.loc[~raw_df.index.isin(df_meta.index)] if len(df_meta) != len(raw_df) else pd.DataFrame()
    raw_df.to_csv(outdir / "tables" / "input_dataset_copy.csv", index=False)
    df_meta.to_csv(outdir / "tables" / "validated_molecules.csv", index=False)
    if not invalid_df.empty:
        invalid_df.to_csv(outdir / "tables" / "invalid_molecules.csv", index=False)

    print(f"Valid molecules: {len(df_meta)} / {len(raw_df)}")
    print("Calculating RDKit physicochemical descriptors...")
    physchem = calc_physchem_descriptors(mols)
    save_table_with_id(df_meta, physchem, outdir / "tables" / "descriptors_physchem_rdkit.csv")

    print("Calculating MACCS fingerprints...")
    maccs = calc_maccs_fingerprints(mols)
    save_table_with_id(df_meta, maccs, outdir / "tables" / "fingerprints_maccs.csv")

    print("Calculating Morgan fingerprints...")
    morgan = calc_morgan_fingerprints(mols, radius=args.morgan_radius, n_bits=args.morgan_bits)
    morgan_name = f"morgan_r{args.morgan_radius}_{args.morgan_bits}"
    save_table_with_id(df_meta, morgan, outdir / "tables" / f"fingerprints_{morgan_name}.csv")

    print("Saving group descriptor summaries...")
    make_group_descriptor_means(df_meta, physchem, outdir)

    print("Running PCA...")
    evr_rows = []
    _, evr = run_pca(physchem.values, df_meta, "physchem", outdir, scale=True)
    evr_rows.append({"representation": "physchem", "PC1_explained_variance": evr[0], "PC2_explained_variance": evr[1], "PC1_PC2_total": sum(evr)})

    _, evr = run_pca(maccs.values, df_meta, "maccs", outdir, scale=False)
    evr_rows.append({"representation": "maccs", "PC1_explained_variance": evr[0], "PC2_explained_variance": evr[1], "PC1_PC2_total": sum(evr)})

    _, evr = run_pca(morgan.values, df_meta, morgan_name, outdir, scale=False)
    evr_rows.append({"representation": morgan_name, "PC1_explained_variance": evr[0], "PC2_explained_variance": evr[1], "PC1_PC2_total": sum(evr)})

    evr_table = pd.DataFrame(evr_rows)
    evr_table.to_csv(outdir / "tables" / "pca_explained_variance.csv", index=False)

    if not args.skip_tsne:
        print("Running t-SNE with selected perplexity values...")
        run_tsne(physchem.values, df_meta, "physchem", outdir, metric="euclidean", perplexities=tsne_perplexities, scale=True, max_iter=args.tsne_iterations)
        run_tsne(maccs.values, df_meta, "maccs", outdir, metric="precomputed_jaccard", perplexities=tsne_perplexities, scale=False, max_iter=args.tsne_iterations)
        run_tsne(morgan.values, df_meta, morgan_name, outdir, metric="precomputed_jaccard", perplexities=tsne_perplexities, scale=False, max_iter=args.tsne_iterations)
    else:
        print("Skipping t-SNE (--skip-tsne).")

    print("Running optional UMAP if available...")
    run_umap_optional(physchem.values, df_meta, "physchem", outdir, metric="euclidean", scale=True)
    run_umap_optional(maccs.values, df_meta, "maccs", outdir, metric="jaccard", scale=False)
    run_umap_optional(morgan.values, df_meta, morgan_name, outdir, metric="jaccard", scale=False)

    print("Running hierarchical clustering...")
    metrics_rows = []
    _, metrics = run_hierarchical_clustering(physchem.values, df_meta, "physchem", outdir, kind="physchem", scale=True)
    metrics_rows.append({"representation": "physchem", **metrics})
    _, metrics = run_hierarchical_clustering(maccs.values, df_meta, "maccs", outdir, kind="fingerprint", scale=False)
    metrics_rows.append({"representation": "maccs", **metrics})
    _, metrics = run_hierarchical_clustering(morgan.values, df_meta, morgan_name, outdir, kind="fingerprint", scale=False)
    metrics_rows.append({"representation": morgan_name, **metrics})

    metrics_table = pd.DataFrame(metrics_rows)
    metrics_table.to_csv(outdir / "tables" / "clustering_separability_metrics.csv", index=False)

    if not args.skip_images:
        print("Drawing 2D molecule grids...")
        draw_molecule_grids(df_meta, mols, outdir)

    sdf_status = None
    if not args.skip_3d:
        print("Generating 3D SDF for Jmol / Avogadro...")
        sdf_status = write_3d_sdf(df_meta, mols, outdir)

    print("Writing summary report...")
    write_summary_report(df_meta, outdir, evr_table, metrics_table, umap_available=umap is not None, sdf_status=sdf_status)

    print("Done.")
    print(f"Results saved to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
