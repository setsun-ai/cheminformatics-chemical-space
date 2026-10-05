# Cheminformatics: exploring chemical space with descriptors, fingerprints and unsupervised learning

Three related cheminformatics projects:
- RDKit descriptors and fingerprints
- dimensionality reduction (PCA, t-SNE, UMAP)
- hierarchical clustering
- comparing how well different molecular representations separate known groups

MSc student project, Gdańsk University of Technology, 2026. Course: *Cheminformatics / chemometrics*.

![Representative molecules](fragrance/results/structures_2d/representative_molecules_grid.png)

## 1. Fragrance and flavour molecules (`fragrance/`)

90 molecules in 3 odour classes (floral, fruity, spicy; 30 each), compared under three representations: physicochemical descriptors, MACCS keys and Morgan fingerprints (r = 2, 2048 bits).

| Representation | PC1 + PC2 variance | ARI (HCA k = 3 vs class) | Silhouette (true classes) |
|---|---|---|---|
| Physchem (RDKit, scaled) | 61.9 % | 0.158 | 0.079 |
| MACCS (Jaccard) | 39.9 % | 0.001 | 0.100 |
| Morgan r2 (Jaccard) | 20.1 % | 0.014 | 0.061 |

Odour class is only weakly encoded by structure. Physicochemical descriptors separate the classes best.

The pipeline also produces:
- t-SNE plots (perplexity 5/15/30)
- UMAP grids
- dendrograms
- a 3D SDF for Jmol

<p>
<img src="fragrance/results/plots_pca/pca_physchem.png" width="49%">
<img src="fragrance/results/plots_hclust/hclust_dendrogram_physchem.png" width="49%">
</p>

```bash
pip install -r requirements.txt
python fragrance/chemoinf_fragrance_analysis.py --input fragrance/data/fragrance_dataset.csv --outdir fragrance/results
```

## 2. "Home pharmacy" dataset (`home_pharmacy/`)

90 active pharmaceutical ingredients in 3 therapeutic groups (pain, respiratory, gastro). The script:
- fetches SMILES/CIDs from **PubChem**
- standardises structures with RDKit (salt stripping, parent fragment, uncharging)
- computes continuous and topological descriptors plus Morgan/MACCS fingerprints
- writes Orange-ready CSVs

The Orange workflow is in `orange/home_pharmacy_projekt.ows`.

```bash
python home_pharmacy/home_pharmacy_cheminfo_project_v5.py --outdir home_pharmacy/data_home_pharmacy_v5
```

## 3. PCA and HCA of the Wine dataset (`reports/PCA_HCA_report.pdf`)

Unsupervised analysis in Orange Data Mining:
- normalisation
- HCA with Euclidean distance (Average / Complete / Ward linkage)
- PCA with explained variance, score plot and loadings

`exploration/` holds earlier descriptor and PubChem scripts.

**Data sources:** all public.
- **Structures:** SMILES/CIDs from [PubChem](https://pubchem.ncbi.nlm.nih.gov/); the fetch log and cache are included.
- **Wine dataset:** the classic UCI Wine dataset, which ships with Orange.

**Tools:** Python (RDKit, scikit-learn, umap-learn, SciPy, matplotlib), PubChem PUG-REST, Orange Data Mining, Jmol.

## AI assistance

The code in this repository was written with the help of AI tools (large language models). Defining the tasks, running the analyses, and checking and interpreting the results were my part of the work.

---

## 🇵🇱 Opis po polsku

Projekty z *Chemoinformatyki*:
1. **Zapachy:** 90 cząsteczek w 3 klasach zapachowych. Porównanie deskryptorów fizykochemicznych, MACCS i Morgan w PCA, t-SNE, UMAP i HCA.
2. **Domowa apteczka:** 90 substancji czynnych pobranych z PubChem, standaryzacja w RDKit, deskryptory i fingerprinty do analizy w Orange.
3. **Raport PCA/HCA** dla zbioru Wine, wykonany w Orange.

Projekt studencki (studia II stopnia), Politechnika Gdańska, 2026.

**Wsparcie AI:** kod w tym repozytorium powstał z pomocą narzędzi AI (dużych modeli językowych). Określenie zadań, uruchamianie analiz oraz sprawdzenie i interpretacja wyników były moją częścią pracy.
