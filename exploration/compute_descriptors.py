"""
compute_descriptors.py
======================
Czyta data/zapachy.csv (output z fetch_and_describe.py) i generuje 3 osobne
pliki CSV, kazdy z tymi samymi 90 wierszami ale innym zestawem deskryptorow:

  data/desc_D1_physchem.csv  - ciagle deskryptory fizykochemiczne (RDKit, ~15 kol.)
  data/desc_D2_morganFP.csv  - Morgan fingerprints (ECFP4, 1024 bity)
  data/desc_D3_maccs.csv     - MACCS keys (166 bitow, interpretowalne podstruktury)

Kazdy plik zaczyna sie od kolumn meta (group, subgroup, name, cid, smiles) oraz
"target variables" (logP, mw, rotb, hbd, heavy_atoms) - w Orange'u oznaczamy je
jako role 'target'/'meta', a kolumny deskryptorow jako 'feature'.

Wymagania: rdkit, ten sam env "chemo".
Uruchomienie: w run.bat zamien linijke python na compute_descriptors.py
albo odpal: python compute_descriptors.py
"""

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, MACCSkeys, rdMolDescriptors

# --- Konfiguracja ---
BASE = Path(__file__).parent
IN_CSV  = BASE / "data" / "zapachy.csv"
OUT_DIR = BASE / "data"
MORGAN_RADIUS = 2     # ECFP4 = radius 2
MORGAN_BITS   = 1024  # 1024 bity to standardowy kompromis - duzo info, ale Orange to udzwignie

# --- Kolumny wspolne (meta + target variables) ---
META_COLS   = ["group", "subgroup", "name", "cid", "smiles"]
TARGET_COLS = ["mw", "logP", "tpsa", "hbd", "hba", "rotb", "heavy_atoms", "rings", "aromatic_rings"]


def physchem_descriptors(mol) -> dict:
    """D1: ciagle deskryptory fizykochemiczne (rozszerzone vs to co bylo w fetch_and_describe)."""
    return {
        # juz mamy w CSV (powtarzamy zeby plik byl samodzielny):
        "MW":             round(Descriptors.MolWt(mol), 3),
        "logP":           round(Crippen.MolLogP(mol), 3),
        "TPSA":           round(Descriptors.TPSA(mol), 3),
        "HBD":            Lipinski.NumHDonors(mol),
        "HBA":            Lipinski.NumHAcceptors(mol),
        "RotB":           Lipinski.NumRotatableBonds(mol),
        "HeavyAtoms":     mol.GetNumHeavyAtoms(),
        "Rings":          rdMolDescriptors.CalcNumRings(mol),
        "AromaticRings":  rdMolDescriptors.CalcNumAromaticRings(mol),
        # dodatkowe, ciekawe:
        "FractionCSP3":   round(rdMolDescriptors.CalcFractionCSP3(mol), 3),  # ile wegli sp3 / wszystkich C
        "NumAliphaticRings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "NumSaturatedRings": rdMolDescriptors.CalcNumSaturatedRings(mol),
        "NumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "MolMR":          round(Crippen.MolMR(mol), 3),     # refrakcja molowa (~rozmiar polaryzowalny)
        "BalabanJ":       round(Descriptors.BalabanJ(mol), 4),    # indeks topologiczny
        "BertzCT":        round(Descriptors.BertzCT(mol), 3),     # zlozonosc topologiczna
        "LabuteASA":      round(Descriptors.LabuteASA(mol), 3),   # powierzchnia dostepna rozpuszczalnikowi
    }


def morgan_fp_bits(mol) -> list[int]:
    """D2: Morgan FP (ECFP4), 1024 bity. Zwraca liste 0/1."""
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=MORGAN_RADIUS, nBits=MORGAN_BITS)
    return list(fp)


def maccs_bits(mol) -> list[int]:
    """D3: MACCS keys, 166 bitow (RDKit zwraca 167 - pomijamy bit 0)."""
    fp = MACCSkeys.GenMACCSKeys(mol)
    bits = list(fp)
    return bits[1:]  # bit 0 nie jest uzywany


def main():
    rows = list(csv.DictReader(IN_CSV.open(encoding="utf-8")))
    print(f"Wczytano {len(rows)} zwiazkow z {IN_CSV.name}")

    # przygotuj molekuly raz
    mols = []
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            raise RuntimeError(f"RDKit nie sparsowal SMILES dla {r['name']}: {r['smiles']}")
        mols.append(m)

    # --- D1: physchem ---
    print("\n[D1] Licze deskryptory fizykochemiczne...")
    physchem = [physchem_descriptors(m) for m in mols]
    pc_cols = list(physchem[0].keys())
    d1_path = OUT_DIR / "desc_D1_physchem.csv"
    with d1_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(META_COLS + pc_cols)
        for r, pc in zip(rows, physchem):
            w.writerow([r[c] for c in META_COLS] + [pc[c] for c in pc_cols])
    print(f"   -> {d1_path.name}: {len(pc_cols)} kolumn ciaglych")

    # --- D2: Morgan FP ---
    print(f"\n[D2] Licze Morgan FP (ECFP4, {MORGAN_BITS} bitow)...")
    fps = [morgan_fp_bits(m) for m in mols]
    fp_cols = [f"FP_{i}" for i in range(MORGAN_BITS)]
    d2_path = OUT_DIR / "desc_D2_morganFP.csv"
    with d2_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(META_COLS + fp_cols)
        for r, fp in zip(rows, fps):
            w.writerow([r[c] for c in META_COLS] + fp)
    # statystyka: ile bitow nieaktywnych w calym zbiorze (mozna potem odrzucic)
    active_per_bit = [sum(fp[i] for fp in fps) for i in range(MORGAN_BITS)]
    n_dead = sum(1 for x in active_per_bit if x == 0)
    n_universal = sum(1 for x in active_per_bit if x == 90)
    print(f"   -> {d2_path.name}: {MORGAN_BITS} bitow")
    print(f"      ({n_dead} bitow zawsze=0 i {n_universal} zawsze=1 - mozesz je usunac w Orange filtrami)")

    # --- D3: MACCS ---
    print("\n[D3] Licze MACCS keys (166 bitow)...")
    maccs = [maccs_bits(m) for m in mols]
    mk_cols = [f"MACCS_{i+1}" for i in range(166)]
    d3_path = OUT_DIR / "desc_D3_maccs.csv"
    with d3_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(META_COLS + mk_cols)
        for r, mk in zip(rows, maccs):
            w.writerow([r[c] for c in META_COLS] + mk)
    active_maccs = [sum(mk[i] for mk in maccs) for i in range(166)]
    n_dead_maccs = sum(1 for x in active_maccs if x == 0)
    print(f"   -> {d3_path.name}: 166 bitow ({n_dead_maccs} zawsze=0)")

    print("\n=================================================")
    print(" Gotowe. Pliki do wgrania do Orange:")
    print(f"   {d1_path}")
    print(f"   {d2_path}")
    print(f"   {d3_path}")
    print("=================================================")


if __name__ == "__main__":
    main()
