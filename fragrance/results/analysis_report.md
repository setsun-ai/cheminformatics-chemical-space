# Chemoinformatics fragrance analysis — automatic report

## Dataset

- Valid molecules: **90**
- Duplicated canonical SMILES: **0**
- Group counts:
  - Floral: 30
  - Fruity: 30
  - Spicy: 30

## Descriptor sets

- Physchem: scaled RDKit physicochemical descriptors for PCA/t-SNE/HCA.
- MACCS: 167-bit structural keys, analysed mainly with Jaccard distance.
- Morgan: radius 2, 2048-bit circular fingerprints, analysed mainly with Jaccard distance.

## PCA explained variance

| representation   |   PC1_explained_variance |   PC2_explained_variance |   PC1_PC2_total |
|:-----------------|-------------------------:|-------------------------:|----------------:|
| physchem         |                 0.39508  |                0.224053  |        0.619133 |
| maccs            |                 0.266121 |                0.133182  |        0.399303 |
| morgan_r2_2048   |                 0.116778 |                0.0843441 |        0.201122 |

## Clustering / separability metrics

| representation   |   ARI_HCA_k3_vs_group |   NMI_HCA_k3_vs_group |   silhouette_true_groups |   silhouette_HCA_k3 |
|:-----------------|----------------------:|----------------------:|-------------------------:|--------------------:|
| physchem         |           0.157894    |             0.166963  |                0.0793555 |           0.213393  |
| maccs            |           0.000534144 |             0.0584587 |                0.099737  |           0.252012  |
| morgan_r2_2048   |           0.0143175   |             0.101628  |                0.0609201 |           0.0820747 |

## UMAP

UMAP plots were generated because `umap-learn` is installed.

## 3D SDF for Jmol

Generated 3D conformers for **90/90** molecules.
Open `structures_3d/fragrance_molecules_3d_for_jmol.sdf` in Jmol/Avogadro for 3D inspection.

## Recommended first plots for presentation

1. `structures_2d/representative_molecules_grid.png`
2. `plots_pca/pca_physchem.png`
3. `plots_pca/pca_maccs.png`
4. `plots_pca/pca_morgan_r2_2048.png`
5. one t-SNE plot for Morgan, preferably perplexity 15 or 30
6. `plots_hclust/hclust_dendrogram_morgan_r2_2048.png`