"""Smoke tests: dataset integrity, a fast run of the fragrance pipeline, offline rebuild of the home-pharmacy data."""
import importlib.util
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = mod
    spec.loader.exec_module(mod)
    return mod


def test_fragrance_dataset_design():
    df = pd.read_csv(ROOT / "fragrance" / "data" / "fragrance_dataset.csv")
    assert df.groupby("group").size().to_dict() == {"Floral": 30, "Fruity": 30, "Spicy": 30}
    canon = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in df["smiles"]]
    assert all(canon) and len(set(canon)) == len(canon)  # all valid, no duplicate structures


def test_fragrance_pipeline_quick_run(tmp_path):
    df = pd.read_csv(ROOT / "fragrance" / "data" / "fragrance_dataset.csv")
    subset = df.groupby("group").head(6)  # 18 molecules, 6 per group
    subset.to_csv(tmp_path / "subset.csv", index=False)
    mod = load(ROOT / "fragrance" / "chemoinf_fragrance_analysis.py")
    mod.main(["--input", str(tmp_path / "subset.csv"), "--outdir", str(tmp_path / "out"),
              "--skip-tsne", "--skip-3d", "--skip-images"])
    tables = tmp_path / "out" / "tables"
    var = pd.read_csv(tables / "pca_explained_variance.csv")
    assert (var["PC1_PC2_total"] <= 1.0 + 1e-9).all()
    metrics = pd.read_csv(tables / "clustering_separability_metrics.csv")
    assert metrics["ARI_HCA_k3_vs_group"].between(-1, 1).all()
    assert (tmp_path / "out" / "analysis_report.md").stat().st_size > 0


def test_home_pharmacy_rebuilds_offline_from_cache(tmp_path, monkeypatch):
    """With the committed PubChem cache the script must not touch the network and must reproduce the CSVs."""
    src = ROOT / "home_pharmacy" / "data_home_pharmacy_v5"
    shutil.copy(src / "pubchem_cache.json", tmp_path / "pubchem_cache.json")
    mod = load(ROOT / "home_pharmacy" / "home_pharmacy_cheminfo_project_v5.py")

    def no_network(*args, **kwargs):
        raise AssertionError("network access attempted")

    monkeypatch.setattr(mod.requests, "get", no_network)
    mod.main(["--outdir", str(tmp_path), "--sleep", "0"])
    new = pd.read_csv(tmp_path / "01_master_descriptors.csv")
    old = pd.read_csv(src / "01_master_descriptors.csv")
    pd.testing.assert_frame_equal(new, old, check_exact=False, rtol=1e-9)


@pytest.mark.parametrize("name", ["02_continuous_for_orange.csv", "05_maccs166_for_orange.csv"])
def test_home_pharmacy_orange_exports_have_90_rows(name):
    df = pd.read_csv(ROOT / "home_pharmacy" / "data_home_pharmacy_v5" / name)
    assert len(df) == 90
