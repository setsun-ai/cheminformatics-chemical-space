"""
fetch_and_describe.py
=====================
Pobiera SMILES z PubChem dla 90 zapachow (3 grupy x 30), waliduje RDKit-em
i liczy deskryptory niezalezne potrzebne jako 'target variables' w Orange'u.

Wyjscie:
    data/zapachy.csv     -- glowny plik do Orange'a
    data/zapachy_log.txt -- log + lista pozycji, ktore wymagaly recznej korekty

Wymagania:
    pip install rdkit requests
    (lub: conda install -c conda-forge rdkit requests)

Uruchomienie:
    python fetch_and_describe.py
"""

import csv
import json
import time
from pathlib import Path
from urllib.parse import quote

import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
SLEEP   = 0.25  # zgodnie z polityka PubChem (max ~5 req/s)

# ---------------------------------------------------------------------------
# 90 zwiazkow w 3 grupach. Format: (group, subgroup, name).
# Nazwy sa dobrane tak, zeby PubChem REST je rozpoznal po polu "name".
# Mosty miedzy grupami sa oznaczone komentarzem [BRIDGE].
# ---------------------------------------------------------------------------

COMPOUNDS = [
    # ============== GRUPA A: TERPENY I TERPENOIDY ==============
    # --- A1: monoterpeny acykliczne ---
    ("A_terpeny", "acyclic_mono", "beta-myrcene"),
    ("A_terpeny", "acyclic_mono", "beta-ocimene"),
    ("A_terpeny", "acyclic_mono", "linalool"),           # [BRIDGE -> octan linalylu w B]
    ("A_terpeny", "acyclic_mono", "geraniol"),           # [BRIDGE -> octan geranylu w B, citral w C]
    ("A_terpeny", "acyclic_mono", "nerol"),
    ("A_terpeny", "acyclic_mono", "citronellol"),        # [BRIDGE -> octan citronellylu w B, citronellal w C]
    # --- A2: monoterpeny monocykliczne ---
    ("A_terpeny", "monocyclic",   "limonene"),
    ("A_terpeny", "monocyclic",   "alpha-terpineol"),
    ("A_terpeny", "monocyclic",   "gamma-terpinene"),
    ("A_terpeny", "monocyclic",   "terpinolene"),
    ("A_terpeny", "monocyclic",   "menthol"),
    ("A_terpeny", "monocyclic",   "carvone"),
    ("A_terpeny", "monocyclic",   "pulegone"),
    ("A_terpeny", "monocyclic",   "p-cymene"),
    # --- A3: monoterpeny bicykliczne ---
    ("A_terpeny", "bicyclic",     "alpha-pinene"),
    ("A_terpeny", "bicyclic",     "beta-pinene"),
    ("A_terpeny", "bicyclic",     "camphor"),
    ("A_terpeny", "bicyclic",     "borneol"),
    ("A_terpeny", "bicyclic",     "sabinene"),
    ("A_terpeny", "bicyclic",     "3-carene"),
    ("A_terpeny", "bicyclic",     "eucalyptol"),         # = 1,8-cineole
    ("A_terpeny", "bicyclic",     "fenchone"),
    # --- A4: seskwiterpeny ---
    ("A_terpeny", "sesqui",       "beta-caryophyllene"),
    ("A_terpeny", "sesqui",       "alpha-humulene"),
    ("A_terpeny", "sesqui",       "farnesol"),
    ("A_terpeny", "sesqui",       "nerolidol"),
    ("A_terpeny", "sesqui",       "alpha-bisabolol"),
    ("A_terpeny", "sesqui",       "cedrol"),
    ("A_terpeny", "sesqui",       "guaiol"),
    ("A_terpeny", "sesqui",       "valencene"),

    # ============== GRUPA B: ESTRY OWOCOWE ==============
    # --- B1: octany krotkie alifatyczne ---
    ("B_estry",   "acetate_short",   "ethyl acetate"),
    ("B_estry",   "acetate_short",   "propyl acetate"),
    ("B_estry",   "acetate_short",   "butyl acetate"),
    ("B_estry",   "acetate_short",   "isoamyl acetate"),
    ("B_estry",   "acetate_short",   "hexyl acetate"),
    ("B_estry",   "acetate_short",   "isobutyl acetate"),
    # --- B2: octany terpenowe/aromatyczne (MOSTY do A) ---
    ("B_estry",   "acetate_terpene_aromatic", "benzyl acetate"),       # [BRIDGE -> benzaldehyd w C]
    ("B_estry",   "acetate_terpene_aromatic", "linalyl acetate"),      # [BRIDGE]
    ("B_estry",   "acetate_terpene_aromatic", "geranyl acetate"),      # [BRIDGE]
    ("B_estry",   "acetate_terpene_aromatic", "citronellyl acetate"),  # [BRIDGE]
    ("B_estry",   "acetate_terpene_aromatic", "phenethyl acetate"),
    ("B_estry",   "acetate_terpene_aromatic", "cinnamyl acetate"),     # [BRIDGE -> cinnamaldehyd w C]
    # --- B3: maslany ---
    ("B_estry",   "butyrate",     "methyl butyrate"),
    ("B_estry",   "butyrate",     "ethyl butyrate"),
    ("B_estry",   "butyrate",     "propyl butyrate"),
    ("B_estry",   "butyrate",     "butyl butyrate"),
    ("B_estry",   "butyrate",     "isoamyl butyrate"),
    ("B_estry",   "butyrate",     "hexyl butyrate"),
    ("B_estry",   "butyrate",     "amyl butyrate"),
    # --- B4: dluzsze lancuchy ---
    ("B_estry",   "long_chain",   "methyl hexanoate"),
    ("B_estry",   "long_chain",   "ethyl hexanoate"),
    ("B_estry",   "long_chain",   "methyl octanoate"),
    ("B_estry",   "long_chain",   "ethyl octanoate"),
    ("B_estry",   "long_chain",   "ethyl decanoate"),
    # --- B5: aromatyczne estry (nie octanowe) ---
    ("B_estry",   "aromatic_ester","methyl salicylate"),    # [BRIDGE -> salicylaldehyd w C]
    ("B_estry",   "aromatic_ester","methyl anthranilate"),
    ("B_estry",   "aromatic_ester","ethyl cinnamate"),      # [BRIDGE -> cinnamaldehyd w C]
    ("B_estry",   "aromatic_ester","benzyl benzoate"),
    # --- B6: inne male estry ---
    ("B_estry",   "other_small",  "ethyl propionate"),
    ("B_estry",   "other_small",  "ethyl formate"),

    # ============== GRUPA C: ALDEHYDY ZAPACHOWE ==============
    # --- C1: alifatyczne nasycone C6-C14 (klasyczne "C-aldehydy") ---
    ("C_aldehydy","aliphatic_sat","hexanal"),
    ("C_aldehydy","aliphatic_sat","heptanal"),
    ("C_aldehydy","aliphatic_sat","octanal"),
    ("C_aldehydy","aliphatic_sat","nonanal"),
    ("C_aldehydy","aliphatic_sat","decanal"),
    ("C_aldehydy","aliphatic_sat","undecanal"),
    ("C_aldehydy","aliphatic_sat","dodecanal"),
    ("C_aldehydy","aliphatic_sat","tridecanal"),
    ("C_aldehydy","aliphatic_sat","tetradecanal"),
    # --- C2: alifatyczne nienasycone "zielone" ---
    ("C_aldehydy","aliphatic_unsat","trans-2-hexenal"),
    ("C_aldehydy","aliphatic_unsat","cis-3-hexenal"),
    ("C_aldehydy","aliphatic_unsat","10-undecenal"),
    ("C_aldehydy","aliphatic_unsat","trans-2-nonenal"),
    # --- C3: aldehydy aromatyczne ---
    ("C_aldehydy","aromatic",     "benzaldehyde"),         # [BRIDGE]
    ("C_aldehydy","aromatic",     "cinnamaldehyde"),       # [BRIDGE]
    ("C_aldehydy","aromatic",     "anisaldehyde"),
    ("C_aldehydy","aromatic",     "vanillin"),
    ("C_aldehydy","aromatic",     "ethylvanillin"),
    ("C_aldehydy","aromatic",     "piperonal"),
    ("C_aldehydy","aromatic",     "cuminaldehyde"),
    ("C_aldehydy","aromatic",     "salicylaldehyde"),      # [BRIDGE]
    ("C_aldehydy","aromatic",     "p-tolualdehyde"),
    ("C_aldehydy","aromatic",     "phenylacetaldehyde"),
    ("C_aldehydy","aromatic",     "helional"),
    # --- C4: terpenowe aldehydy (MOSTY do A) ---
    ("C_aldehydy","terpene_ald",  "citral"),               # [BRIDGE]
    ("C_aldehydy","terpene_ald",  "citronellal"),          # [BRIDGE]
    ("C_aldehydy","terpene_ald",  "perillaldehyde"),
    ("C_aldehydy","terpene_ald",  "safranal"),
    ("C_aldehydy","terpene_ald",  "myrtenal"),
    ("C_aldehydy","terpene_ald",  "hydroxycitronellal"),
]


# ---------------------------------------------------------------------------
# PubChem API
# ---------------------------------------------------------------------------

def get_cid(name: str) -> int | None:
    """Zwraca pierwszy CID dla nazwy lub None."""
    url = f"{PUBCHEM}/compound/name/{quote(name)}/cids/JSON"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            return None
        cids = r.json().get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except Exception as e:
        print(f"  ! get_cid({name}): {e}")
        return None


def get_smiles(cid: int) -> tuple[str | None, str | None]:
    """Zwraca (canonical_smiles, isomeric_smiles) dla CID."""
    # PubChem zmienial nazwy property: probujemy najpierw 'SMILES', potem stare nazwy
    for prop_str in ["IsomericSMILES,CanonicalSMILES", "SMILES,ConnectivitySMILES", "SMILES"]:
        url = f"{PUBCHEM}/compound/cid/{cid}/property/{prop_str}/JSON"
        try:
            r = requests.get(url, timeout=15)
            if r.status_code != 200:
                continue
            props = r.json().get("PropertyTable", {}).get("Properties", [])
            if not props:
                continue
            p = props[0]
            iso = p.get("IsomericSMILES") or p.get("SMILES")
            can = p.get("CanonicalSMILES") or p.get("ConnectivitySMILES") or iso
            if iso or can:
                return can, iso
        except Exception as e:
            print(f"  ! get_smiles(CID={cid}): {e}")
    return None, None


# ---------------------------------------------------------------------------
# RDKit: walidacja + deskryptory
# ---------------------------------------------------------------------------

def describe(smiles: str) -> dict | None:
    """Liczy deskryptory niezalezne dla SMILES (po walidacji RDKit)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "rdkit_smiles":   Chem.MolToSmiles(mol),
        "mw":             round(Descriptors.MolWt(mol), 3),
        "logP":           round(Crippen.MolLogP(mol), 3),
        "tpsa":           round(Descriptors.TPSA(mol), 3),
        "hbd":            Lipinski.NumHDonors(mol),
        "hba":            Lipinski.NumHAcceptors(mol),
        "rotb":           Lipinski.NumRotatableBonds(mol),
        "heavy_atoms":    mol.GetNumHeavyAtoms(),
        "rings":          rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(exist_ok=True)

    rows = []
    log_lines = []

    print(f"Pobieram {len(COMPOUNDS)} zwiazkow z PubChem...\n")

    for i, (group, subgroup, name) in enumerate(COMPOUNDS, 1):
        print(f"[{i:>2}/{len(COMPOUNDS)}] {group:>11s} / {subgroup:<25s} {name}")
        cid = get_cid(name)
        time.sleep(SLEEP)
        if cid is None:
            msg = f"  -> BRAK CID dla '{name}'"
            print(msg)
            log_lines.append(msg)
            rows.append({"group": group, "subgroup": subgroup, "name": name,
                         "cid": "", "smiles": "", "isomeric_smiles": "",
                         "mw": "", "logP": "", "tpsa": "",
                         "hbd": "", "hba": "", "rotb": "",
                         "heavy_atoms": "", "rings": "", "aromatic_rings": "",
                         "status": "no_cid"})
            continue

        can, iso = get_smiles(cid)
        time.sleep(SLEEP)
        if not (can or iso):
            msg = f"  -> BRAK SMILES dla CID={cid} ({name})"
            print(msg)
            log_lines.append(msg)
            rows.append({"group": group, "subgroup": subgroup, "name": name,
                         "cid": cid, "smiles": "", "isomeric_smiles": "",
                         "mw": "", "logP": "", "tpsa": "",
                         "hbd": "", "hba": "", "rotb": "",
                         "heavy_atoms": "", "rings": "", "aromatic_rings": "",
                         "status": "no_smiles"})
            continue

        smiles_for_rdkit = iso or can
        d = describe(smiles_for_rdkit)
        if d is None:
            msg = f"  -> RDKit nie sparsowal SMILES: '{smiles_for_rdkit}' ({name}, CID={cid})"
            print(msg)
            log_lines.append(msg)
            rows.append({"group": group, "subgroup": subgroup, "name": name,
                         "cid": cid, "smiles": can or "", "isomeric_smiles": iso or "",
                         "mw": "", "logP": "", "tpsa": "",
                         "hbd": "", "hba": "", "rotb": "",
                         "heavy_atoms": "", "rings": "", "aromatic_rings": "",
                         "status": "rdkit_fail"})
            continue

        rows.append({
            "group":          group,
            "subgroup":       subgroup,
            "name":           name,
            "cid":            cid,
            "smiles":         d["rdkit_smiles"],   # kanoniczny SMILES z RDKit (ujednolicony)
            "isomeric_smiles": iso or "",
            "mw":             d["mw"],
            "logP":           d["logP"],
            "tpsa":           d["tpsa"],
            "hbd":            d["hbd"],
            "hba":            d["hba"],
            "rotb":           d["rotb"],
            "heavy_atoms":    d["heavy_atoms"],
            "rings":          d["rings"],
            "aromatic_rings": d["aromatic_rings"],
            "status":         "ok",
        })

    # --- CSV ---
    csv_path = out_dir / "zapachy.csv"
    fieldnames = ["group", "subgroup", "name", "cid",
                  "smiles", "isomeric_smiles",
                  "mw", "logP", "tpsa", "hbd", "hba", "rotb",
                  "heavy_atoms", "rings", "aromatic_rings", "status"]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\n>>> Zapisano: {csv_path}")

    # --- log ---
    log_path = out_dir / "zapachy_log.txt"
    ok    = sum(1 for r in rows if r["status"] == "ok")
    fails = [r for r in rows if r["status"] != "ok"]
    log_text = (
        f"OK: {ok}/{len(rows)}\n"
        f"Bledy: {len(fails)}\n\n"
        + "\n".join(log_lines)
        + ("\n\nDo recznej korekty:\n" + "\n".join(f"  - {r['name']} ({r['status']})" for r in fails) if fails else "")
    )
    log_path.write_text(log_text, encoding="utf-8")
    print(f">>> Log:       {log_path}")
    print(f"\nOK: {ok}/{len(rows)}")
    if fails:
        print("Do recznej korekty:")
        for r in fails:
            print(f"  - {r['name']:<25s}  status={r['status']}")


if __name__ == "__main__":
    main()
