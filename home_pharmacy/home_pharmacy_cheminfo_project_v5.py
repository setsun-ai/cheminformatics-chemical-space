#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chemoinformatyka — projekt: "domowa apteczka"

Cel skryptu:
1) tworzy listę 90 substancji czynnych w 3 jednowyrazowych grupach:
   Bol, Oddech, Gastro;
2) pobiera SMILES i CID z PubChem;
3) standaryzuje cząsteczki RDKit-em: cleanup, largest/fragment parent, uncharge;
4) liczy deskryptory ciągłe i topologiczne;
5) generuje fingerprinty Morgan/ECFP4 oraz MACCS;
6) zapisuje pliki CSV gotowe do Orange.

Instalacja, najlepiej przez conda:
    conda create -n chemo python=3.11 -y
    conda activate chemo
    conda install -c conda-forge rdkit pandas requests -y

Uruchomienie:
    python home_pharmacy_cheminfo_project_v5.py --outdir data_home_pharmacy_v5

Uwaga:
- Skrypt korzysta z PubChem, więc potrzebuje internetu przy pierwszym uruchomieniu.
- PubChem czasem zwraca sole. Skrypt próbuje sprowadzić je do parent/active moiety.
- Po uruchomieniu sprawdź plik 00_fetch_status.csv i ewentualne rekordy FAILED.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, GraphDescriptors, Lipinski, MACCSkeys, rdMolDescriptors

try:
    from rdkit.Chem import rdFingerprintGenerator
except Exception:  # older RDKit fallback
    rdFingerprintGenerator = None
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass(frozen=True)
class CompoundSeed:
    name: str
    group: str          # one-word group: Bol / Oddech / Gastro
    subgroup: str       # more detailed mechanistic/use subgroup
    query: Optional[str] = None  # optional PubChem query override


COMPOUNDS: List[CompoundSeed] = [
    # ===== Grupa 1: Bol =====
    CompoundSeed("acetaminophen", "Bol", "anilide"),
    CompoundSeed("ibuprofen", "Bol", "nsaid_propionic"),
    CompoundSeed("naproxen", "Bol", "nsaid_propionic"),
    CompoundSeed("ketoprofen", "Bol", "nsaid_propionic"),
    CompoundSeed("dexketoprofen", "Bol", "nsaid_propionic"),
    CompoundSeed("flurbiprofen", "Bol", "nsaid_propionic"),
    CompoundSeed("fenoprofen", "Bol", "nsaid_propionic"),
    CompoundSeed("tiaprofenic acid", "Bol", "nsaid_propionic"),
    CompoundSeed("aspirin", "Bol", "nsaid_salicylate"),
    CompoundSeed("salicylic acid", "Bol", "nsaid_salicylate"),
    CompoundSeed("diclofenac", "Bol", "nsaid_acetic"),
    CompoundSeed("aceclofenac", "Bol", "nsaid_acetic"),
    CompoundSeed("indomethacin", "Bol", "nsaid_acetic"),
    CompoundSeed("ketorolac", "Bol", "nsaid_acetic"),
    CompoundSeed("mefenamic acid", "Bol", "nsaid_fenamate"),
    CompoundSeed("tolfenamic acid", "Bol", "nsaid_fenamate"),
    CompoundSeed("meloxicam", "Bol", "nsaid_oxicam"),
    CompoundSeed("piroxicam", "Bol", "nsaid_oxicam"),
    CompoundSeed("lornoxicam", "Bol", "nsaid_oxicam"),
    CompoundSeed("nimesulide", "Bol", "nsaid_sulfonanilide"),
    CompoundSeed("celecoxib", "Bol", "nsaid_coxib"),
    CompoundSeed("etoricoxib", "Bol", "nsaid_coxib"),
    CompoundSeed("metamizole", "Bol", "pyrazolone"),
    CompoundSeed("propyphenazone", "Bol", "pyrazolone"),
    CompoundSeed("phenazone", "Bol", "pyrazolone", query="antipyrine"),
    CompoundSeed("caffeine", "Bol", "adjuvant"),
    CompoundSeed("codeine", "Bol", "opioid_like"),
    CompoundSeed("tramadol", "Bol", "opioid_like"),
    CompoundSeed("benzydamine", "Bol", "local_antiinflammatory"),
    CompoundSeed("benzocaine", "Bol", "local_anesthetic"),

    # ===== Grupa 2: Oddech =====
    CompoundSeed("cetirizine", "Oddech", "antihistamine"),
    CompoundSeed("levocetirizine", "Oddech", "antihistamine"),
    CompoundSeed("loratadine", "Oddech", "antihistamine"),
    CompoundSeed("desloratadine", "Oddech", "antihistamine"),
    CompoundSeed("fexofenadine", "Oddech", "antihistamine"),
    CompoundSeed("bilastine", "Oddech", "antihistamine"),
    CompoundSeed("ebastine", "Oddech", "antihistamine"),
    CompoundSeed("rupatadine", "Oddech", "antihistamine"),
    CompoundSeed("diphenhydramine", "Oddech", "antihistamine"),
    CompoundSeed("chlorpheniramine", "Oddech", "antihistamine"),
    CompoundSeed("clemastine", "Oddech", "antihistamine"),
    CompoundSeed("hydroxyzine", "Oddech", "antihistamine"),
    CompoundSeed("phenylephrine", "Oddech", "decongestant"),
    CompoundSeed("pseudoephedrine", "Oddech", "decongestant"),
    CompoundSeed("ephedrine", "Oddech", "decongestant"),
    CompoundSeed("xylometazoline", "Oddech", "decongestant"),
    CompoundSeed("oxymetazoline", "Oddech", "decongestant"),
    CompoundSeed("naphazoline", "Oddech", "decongestant"),
    CompoundSeed("dextromethorphan", "Oddech", "antitussive"),
    CompoundSeed("butamirate", "Oddech", "antitussive", query="butamirate citrate"),
    CompoundSeed("noscapine", "Oddech", "antitussive"),
    CompoundSeed("ambroxol", "Oddech", "mucolytic"),
    CompoundSeed("bromhexine", "Oddech", "mucolytic"),
    CompoundSeed("guaifenesin", "Oddech", "mucolytic"),
    CompoundSeed("acetylcysteine", "Oddech", "mucolytic"),
    CompoundSeed("carbocisteine", "Oddech", "mucolytic", query="carbocysteine"),
    CompoundSeed("erdosteine", "Oddech", "mucolytic"),
    CompoundSeed("salbutamol", "Oddech", "bronchodilator", query="albuterol"),
    CompoundSeed("montelukast", "Oddech", "leukotriene_antagonist"),
    CompoundSeed("azelastine", "Oddech", "nasal_antihistamine"),

    # ===== Grupa 3: Gastro =====
    CompoundSeed("omeprazole", "Gastro", "ppi"),
    CompoundSeed("esomeprazole", "Gastro", "ppi"),
    CompoundSeed("pantoprazole", "Gastro", "ppi"),
    CompoundSeed("lansoprazole", "Gastro", "ppi"),
    CompoundSeed("rabeprazole", "Gastro", "ppi"),
    CompoundSeed("famotidine", "Gastro", "h2_blocker"),
    CompoundSeed("cimetidine", "Gastro", "h2_blocker"),
    CompoundSeed("ranitidine", "Gastro", "h2_blocker"),
    CompoundSeed("nizatidine", "Gastro", "h2_blocker"),
    CompoundSeed("lafutidine", "Gastro", "h2_blocker"),
    CompoundSeed("loperamide", "Gastro", "antidiarrheal"),
    CompoundSeed("diphenoxylate", "Gastro", "antidiarrheal"),
    CompoundSeed("racecadotril", "Gastro", "antidiarrheal"),
    CompoundSeed("metoclopramide", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("domperidone", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("ondansetron", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("granisetron", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("itopride", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("mosapride", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("prucalopride", "Gastro", "antiemetic_prokinetic"),
    CompoundSeed("hyoscine", "Gastro", "antispasmodic", query="scopolamine"),
    CompoundSeed("drotaverine", "Gastro", "antispasmodic"),
    CompoundSeed("mebeverine", "Gastro", "antispasmodic"),
    CompoundSeed("trimebutine", "Gastro", "antispasmodic"),
    CompoundSeed("alverine", "Gastro", "antispasmodic"),
    CompoundSeed("otilonium", "Gastro", "antispasmodic", query="otilonium bromide"),
    CompoundSeed("pinaverium", "Gastro", "antispasmodic", query="pinaverium bromide"),
    CompoundSeed("dicyclomine", "Gastro", "antispasmodic"),
    CompoundSeed("bisacodyl", "Gastro", "laxative"),
    CompoundSeed("picosulfate", "Gastro", "laxative", query="sodium picosulfate"),
]

METADATA_COLS = ["compound_name", "group", "subgroup", "canonical_smiles"]
CONTINUOUS_COLS = [
    "MolWt", "ExactMolWt", "HeavyAtomCount", "MolLogP", "MolMR", "TPSA",
    "HBD", "HBA", "RotBonds", "RingCount", "AromaticRingCount",
    "AliphaticRingCount", "FractionCSP3", "FormalCharge", "NumHeteroatoms",
    "NumAmideBonds", "NumSaturatedRings", "NumAromaticHeterocycles",
    "NumAromaticCarbocycles",
]
TOPOLOGICAL_COLS = [
    "BalabanJ", "BertzCT", "HallKierAlpha", "Ipc", "Kappa1", "Kappa2", "Kappa3",
    "Chi0", "Chi1", "Chi0n", "Chi1n", "Chi0v", "Chi1v",
]
MORGAN_COLS = [f"morgan1024_{i}" for i in range(1024)]
MACCS_COLS = [f"maccs_{i}" for i in range(1, 167)]


# PubChem PUG-REST endpoints.
# PubChem has changed/varied SMILES property names over time. Some installations/API
# responses return CanonicalSMILES/IsomericSMILES, others return
# ConnectivitySMILES/SMILES. Therefore we try multiple property sets and then fall
# back to SDF-by-CID, which RDKit can parse directly.
PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    "{name}/property/{properties}/JSON"
)
PUBCHEM_CIDS_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    "{name}/cids/JSON"
)
PUBCHEM_SDF_BY_CID_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    "{cid}/SDF?record_type=2d"
)

PROPERTY_SETS = [
    "CanonicalSMILES,IsomericSMILES,InChIKey",
    "SMILES,ConnectivitySMILES,InChIKey",
]
SMILES_KEYS_PRIORITY = [
    "IsomericSMILES",
    "CanonicalSMILES",
    "SMILES",
    "ConnectivitySMILES",
]


def http_get(url: str, timeout: int = 30, retries: int = 3) -> requests.Response:
    """Small retry wrapper for PubChem calls.

    It avoids crashing on a transient network hiccup. Final failure is still reported
    in 08_failed_or_problematic.csv rather than stopping the whole run.
    """
    headers = {"User-Agent": "PG-cheminformatics-student-project/1.0"}
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=timeout, headers=headers)
            response.raise_for_status()
            return response
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(0.75 * attempt)
    raise last_error


def _extract_smiles_from_props(props: dict) -> Tuple[Optional[str], Optional[str]]:
    """Return best available SMILES and the property key used."""
    for key in SMILES_KEYS_PRIORITY:
        value = props.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip(), key
    # Defensive fallback: PubChem may introduce another key containing 'SMILES'.
    for key, value in props.items():
        if "smiles" in str(key).lower() and isinstance(value, str) and value.strip():
            return value.strip(), str(key)
    return None, None


def _fetch_first_cid(query: str) -> Optional[int]:
    url = PUBCHEM_CIDS_URL.format(name=quote(query))
    r = http_get(url, timeout=30)
    data = r.json()
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    return int(cids[0])


def _fetch_smiles_from_sdf(cid: int) -> Tuple[Optional[str], Optional[str]]:
    """Fallback: fetch a PubChem 2D SDF record and convert it to SMILES with RDKit."""
    url = PUBCHEM_SDF_BY_CID_URL.format(cid=cid)
    r = http_get(url, timeout=45)
    mol = Chem.MolFromMolBlock(r.text, sanitize=True, removeHs=False)
    if mol is None:
        return None, "SDF parse failed"
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    return smiles, None


def fetch_pubchem_record(query: str, cache: Dict[str, dict], sleep_s: float = 0.20) -> dict:
    """Fetch one PubChem record by compound name, with JSON cache and robust SMILES fallback."""
    key = query.strip().lower()

    # Reuse cache only if it contains a usable SMILES. Older bad-cache records from
    # previous script versions may have status OK but empty SMILES, so they must be retried.
    cached = cache.get(key)
    if cached and cached.get("status") == "OK":
        cached_smiles = cached.get("pubchem_isomeric_smiles") or cached.get("pubchem_canonical_smiles")
        if isinstance(cached_smiles, str) and cached_smiles.strip():
            return cached

    errors = []
    cid = None
    inchikey = None
    smiles = None
    smiles_source = None

    # 1) Try property endpoints with both old and current SMILES property names.
    for prop_set in PROPERTY_SETS:
        url = PUBCHEM_PROPERTY_URL.format(name=quote(query), properties=prop_set)
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            data = r.json()
            props = data["PropertyTable"]["Properties"][0]
            cid = props.get("CID") or cid
            inchikey = props.get("InChIKey") or inchikey
            smiles, smiles_source = _extract_smiles_from_props(props)
            if smiles:
                break
            errors.append(f"No SMILES in property response for {prop_set}; keys={list(props.keys())}")
        except Exception as e:
            errors.append(f"Property endpoint failed for {prop_set}: {repr(e)}")

    # 2) If needed, get CID explicitly.
    if cid is None:
        try:
            cid = _fetch_first_cid(query)
        except Exception as e:
            errors.append(f"CID endpoint failed: {repr(e)}")

    # 3) Final fallback: fetch SDF by CID and convert to SMILES using RDKit.
    if not smiles and cid is not None:
        try:
            sdf_smiles, sdf_error = _fetch_smiles_from_sdf(int(cid))
            if sdf_smiles:
                smiles = sdf_smiles
                smiles_source = "SDF_by_CID"
            elif sdf_error:
                errors.append(sdf_error)
        except Exception as e:
            errors.append(f"SDF endpoint failed: {repr(e)}")

    if smiles:
        record = {
            "status": "OK",
            "query": query,
            "cid": cid,
            # Keep both fields populated to stay compatible with the rest of the script.
            "pubchem_canonical_smiles": smiles,
            "pubchem_isomeric_smiles": smiles,
            "pubchem_inchikey": inchikey,
            "smiles_source": smiles_source or "unknown",
            "error": " | ".join(errors),
        }
    else:
        record = {
            "status": "FAILED",
            "query": query,
            "cid": cid,
            "pubchem_canonical_smiles": None,
            "pubchem_isomeric_smiles": None,
            "pubchem_inchikey": inchikey,
            "smiles_source": None,
            "error": " | ".join(errors) or "No SMILES obtained",
        }

    cache[key] = record
    time.sleep(sleep_s)
    return record


def standardize_smiles(smiles: str) -> Tuple[Optional[Chem.Mol], Optional[str], str]:
    """Return standardized RDKit mol, canonical SMILES, and status text."""
    if not smiles or not isinstance(smiles, str):
        return None, None, "empty_smiles"

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None, "rdkit_parse_failed"

        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)  # removes counterions/small fragments
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        Chem.SanitizeMol(mol)

        can = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return mol, can, "OK"
    except Exception as e:
        return None, None, f"standardization_failed: {repr(e)}"


def safe_float(func, mol: Chem.Mol):
    try:
        if func is None:
            return None
        value = func(mol)
        if value is None:
            return None
        return float(value)
    except Exception:
        return None


def safe_int(func, mol: Chem.Mol):
    try:
        if func is None:
            return None
        value = func(mol)
        if value is None:
            return None
        return int(value)
    except Exception:
        return None


def compute_descriptors(mol: Chem.Mol) -> dict:
    """Continuous + topology-like descriptors useful for Orange/PCA/HCA."""
    return {
        # continuous physicochemical descriptors
        "MolWt": safe_float(Descriptors.MolWt, mol),
        "ExactMolWt": safe_float(Descriptors.ExactMolWt, mol),
        "HeavyAtomCount": safe_int(Descriptors.HeavyAtomCount, mol),
        "MolLogP": safe_float(Crippen.MolLogP, mol),
        "MolMR": safe_float(Crippen.MolMR, mol),
        "TPSA": safe_float(rdMolDescriptors.CalcTPSA, mol),
        "HBD": safe_int(Lipinski.NumHDonors, mol),
        "HBA": safe_int(Lipinski.NumHAcceptors, mol),
        "RotBonds": safe_int(Lipinski.NumRotatableBonds, mol),
        "RingCount": safe_int(rdMolDescriptors.CalcNumRings, mol),
        "AromaticRingCount": safe_int(rdMolDescriptors.CalcNumAromaticRings, mol),
        "AliphaticRingCount": safe_int(rdMolDescriptors.CalcNumAliphaticRings, mol),
        "FractionCSP3": safe_float(rdMolDescriptors.CalcFractionCSP3, mol),
        "FormalCharge": safe_int(Chem.GetFormalCharge, mol),
        "NumHeteroatoms": safe_int(rdMolDescriptors.CalcNumHeteroatoms, mol),
        "NumAmideBonds": safe_int(getattr(rdMolDescriptors, "CalcNumAmideBonds", None), mol),
        "NumSaturatedRings": safe_int(rdMolDescriptors.CalcNumSaturatedRings, mol),
        "NumAromaticHeterocycles": safe_int(rdMolDescriptors.CalcNumAromaticHeterocycles, mol),
        "NumAromaticCarbocycles": safe_int(rdMolDescriptors.CalcNumAromaticCarbocycles, mol),
        # graph/topological descriptors
        "BalabanJ": safe_float(GraphDescriptors.BalabanJ, mol),
        "BertzCT": safe_float(GraphDescriptors.BertzCT, mol),
        "HallKierAlpha": safe_float(GraphDescriptors.HallKierAlpha, mol),
        "Ipc": safe_float(GraphDescriptors.Ipc, mol),
        "Kappa1": safe_float(GraphDescriptors.Kappa1, mol),
        "Kappa2": safe_float(GraphDescriptors.Kappa2, mol),
        "Kappa3": safe_float(GraphDescriptors.Kappa3, mol),
        "Chi0": safe_float(GraphDescriptors.Chi0, mol),
        "Chi1": safe_float(GraphDescriptors.Chi1, mol),
        "Chi0n": safe_float(GraphDescriptors.Chi0n, mol),
        "Chi1n": safe_float(GraphDescriptors.Chi1n, mol),
        "Chi0v": safe_float(GraphDescriptors.Chi0v, mol),
        "Chi1v": safe_float(GraphDescriptors.Chi1v, mol),
    }


def morgan_bits(mol: Chem.Mol, radius: int = 2, n_bits: int = 1024) -> Dict[str, int]:
    """Morgan/ECFP-like bit vector as dict.

    Prefer the modern RDKit MorganGenerator API; fall back to the older
    GetMorganFingerprintAsBitVect API. Uses ToBitString only, so it does not
    depend on NumPy array conversion.
    """
    if rdFingerprintGenerator is not None:
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        fp = generator.GetFingerprint(mol)
    else:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
    bitstring = fp.ToBitString()
    if len(bitstring) != n_bits:
        raise ValueError(f"Unexpected Morgan fingerprint length: {len(bitstring)} != {n_bits}")
    return {f"morgan{n_bits}_{i}": int(bitstring[i]) for i in range(n_bits)}

def maccs_bits(mol: Chem.Mol) -> Dict[str, int]:
    """MACCS keys. RDKit returns 167 bits; bit 0 is not used, so keep 1..166."""
    fp = MACCSkeys.GenMACCSKeys(mol)
    bitstring = fp.ToBitString()
    return {f"maccs_{i}": int(bitstring[i]) for i in range(1, 167)}


def load_cache(path: Path) -> Dict[str, dict]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {}


def save_cache(path: Path, cache: Dict[str, dict]) -> None:
    path.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(Path(__file__).resolve().parent / "data_home_pharmacy_v5"),
                        help="Output directory (also holds pubchem_cache.json)")
    parser.add_argument("--sleep", type=float, default=0.20, help="Delay between PubChem requests")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cache_path = outdir / "pubchem_cache.json"
    cache = load_cache(cache_path)

    rows_master = []
    rows_continuous = []
    rows_topological = []
    rows_morgan = []
    rows_maccs = []
    fetch_status = []

    seen_keys = {}

    for i, seed in enumerate(COMPOUNDS, start=1):
        query = seed.query or seed.name
        print(f"[{i:02d}/{len(COMPOUNDS)}] {seed.group:7s} | {seed.name}  <- PubChem: {query}")

        rec = fetch_pubchem_record(query, cache=cache, sleep_s=args.sleep)
        save_cache(cache_path, cache)

        smiles_raw = rec.get("pubchem_isomeric_smiles") or rec.get("pubchem_canonical_smiles")
        mol, canonical_smiles, std_status = standardize_smiles(smiles_raw)

        base = {
            "compound_name": seed.name,
            "group": seed.group,
            "subgroup": seed.subgroup,
            "pubchem_query": query,
            "pubchem_cid": rec.get("cid"),
            "pubchem_inchikey": rec.get("pubchem_inchikey"),
            "pubchem_smiles": smiles_raw,
            "canonical_smiles": canonical_smiles,
            "standardization_status": std_status,
            "fetch_status": rec.get("status"),
            "fetch_error": rec.get("error", ""),
        }

        fetch_status.append(base.copy())

        if mol is None:
            rows_master.append(base)
            continue

        try:
            std_inchikey = Chem.MolToInchiKey(mol)
        except Exception:
            std_inchikey = canonical_smiles
        base["standardized_inchikey"] = std_inchikey

        if std_inchikey in seen_keys:
            base["duplicate_of"] = seen_keys[std_inchikey]
        else:
            base["duplicate_of"] = ""
            seen_keys[std_inchikey] = seed.name

        desc = compute_descriptors(mol)
        master = {**base, **desc}
        rows_master.append(master)

        metadata = {
            "compound_name": seed.name,
            "group": seed.group,
            "subgroup": seed.subgroup,
            "canonical_smiles": canonical_smiles,
        }

        try:
            rows_continuous.append({**metadata, **{k: desc.get(k) for k in CONTINUOUS_COLS}})
            rows_topological.append({**metadata, **{k: desc.get(k) for k in TOPOLOGICAL_COLS}})
            rows_morgan.append({**metadata, **morgan_bits(mol, radius=2, n_bits=1024)})
            rows_maccs.append({**metadata, **maccs_bits(mol)})
        except Exception as e:
            # Do not crash the whole project because one fingerprint/descriptor failed.
            base["standardization_status"] = f"descriptor_or_fingerprint_failed: {repr(e)}"
            rows_master[-1].update(base)
            continue

    df_master = pd.DataFrame(rows_master)
    df_status = pd.DataFrame(fetch_status)
    df_cont = pd.DataFrame(rows_continuous, columns=METADATA_COLS + CONTINUOUS_COLS)
    df_topo = pd.DataFrame(rows_topological, columns=METADATA_COLS + TOPOLOGICAL_COLS)
    df_morgan = pd.DataFrame(rows_morgan, columns=METADATA_COLS + MORGAN_COLS)
    df_maccs = pd.DataFrame(rows_maccs, columns=METADATA_COLS + MACCS_COLS)

    # Useful QC tables
    group_counts = df_master.groupby("group", dropna=False)["compound_name"].count().reset_index(name="n")
    duplicate_rows = df_master[df_master.get("duplicate_of", "").fillna("") != ""] if "duplicate_of" in df_master.columns else pd.DataFrame()
    failed_rows = df_master[(df_master["fetch_status"] != "OK") | (df_master["standardization_status"] != "OK")]

    # Save outputs
    df_status.to_csv(outdir / "00_fetch_status.csv", index=False)
    df_master.to_csv(outdir / "01_master_descriptors.csv", index=False)
    df_cont.to_csv(outdir / "02_continuous_for_orange.csv", index=False)
    df_topo.to_csv(outdir / "03_topological_for_orange.csv", index=False)
    df_morgan.to_csv(outdir / "04_morgan1024_for_orange.csv", index=False)
    df_maccs.to_csv(outdir / "05_maccs166_for_orange.csv", index=False)
    group_counts.to_csv(outdir / "06_group_counts.csv", index=False)
    duplicate_rows.to_csv(outdir / "07_possible_duplicates.csv", index=False)
    failed_rows.to_csv(outdir / "08_failed_or_problematic.csv", index=False)

    print("\nDone.")
    print(f"Output directory: {outdir.resolve()}")
    print("\nGroup counts:")
    print(group_counts.to_string(index=False))
    print(f"\nUsable standardized molecules: {len(df_cont)} / {len(COMPOUNDS)}")
    if len(df_cont) != len(COMPOUNDS):
        print("WARNING: Not all molecules produced descriptor/fingerprint rows; inspect 08_failed_or_problematic.csv")

    if not failed_rows.empty:
        print("\nWARNING: Some compounds failed or had standardization problems.")
        print(f"Check: {outdir / '08_failed_or_problematic.csv'}")

    if not duplicate_rows.empty:
        print("\nWARNING: Possible duplicates after standardization.")
        print(f"Check: {outdir / '07_possible_duplicates.csv'}")

    print("\nRecommended Orange files:")
    print("- 02_continuous_for_orange.csv  -> PCA/HCA/t-SNE/UMAP with global descriptors")
    print("- 03_topological_for_orange.csv -> topology/scaffold complexity")
    print("- 04_morgan1024_for_orange.csv  -> structural fingerprints")
    print("- 05_maccs166_for_orange.csv    -> interpretable fingerprint-like keys")
    print("\nIn Orange: set 'group' as class/target for coloring; use MolLogP, TPSA, MolWt, HBA/HBD")
    print("as continuous variables for interpretation of trends.")


if __name__ == "__main__":
    main()
