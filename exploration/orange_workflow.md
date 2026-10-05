# Workflow w Orange — krok po kroku

Workflow oparty o metody z wykładu: **PCA + HCA z dendrogramem** (kryterium Sneatha 2/3D, jak na slajdach). k-Means i t-SNE są opcjonalne (bonus „dla ambitnych").

Wszystkie kroki powtarzasz **dla każdego z trzech CSV-ów** (deskryptory ciągłe / Morgan / ChemBERTa) — to jest serce projektu: porównanie reprezentacji.

## 0. Co dostajesz ze skryptu

Folder `output/` zawiera między innymi:

```
combined_continuous.csv   ← deskryptory RDKit 2D (ciągłe) - START STĄD
combined_morgan.csv       ← Morgan ECFP4 (binarne, 1024 bity)
combined_maccs.csv        ← MACCS (binarne, 166 bitów) - opcjonalnie
combined_chemberta.csv    ← embeddingi ChemBERTa (384 ciągłe wymiary)
mol_images/<name>.png     ← obrazki 2D wszystkich cząsteczek
```

Każdy `combined_*.csv` ma tę samą strukturę meta-kolumn na początku:
`name, group, subgroup, MW, logP, TPSA, HBD, HBA, RotB, HeavyAtoms, AromaticRings`
a po nich kolumny z deskryptorami danego typu.

## 1. Wczytaj dane: widget **File**

Drag-and-drop ikonki "File" na canvas, otwórz dwukrotnym kliknięciem.

Wybierz np. `combined_continuous.csv`. W tabeli ról kolumn ustaw:

| Kolumna | Rola |
|---|---|
| `name` | **meta** (etykieta, nie wchodzi do analizy) |
| `group` | **target** (kategoryczna — klasy do kolorowania) |
| `subgroup` | **meta** (do kolorowania jako alternatywa) |
| `MW`, `logP`, `TPSA`, `HBD`, `HBA`, `RotB`, `HeavyAtoms`, `AromaticRings` | **meta** (zostają do interpretacji) |
| reszta (`ExactMW`, `MolMR`, `LabuteASA`, ...) | **feature** |

Po `Apply` powinno być ~28 features dla `combined_continuous.csv`, 1024 dla `combined_morgan.csv`, 384 dla `combined_chemberta.csv`.

> **Dlaczego targety jako meta a nie feature?** Bo nie chcemy, żeby logP/MW „przeciekały" do PCA jako jedne ze zmiennych — to byłoby cyrkulacyjne. Trzymamy je jako meta, żeby później kolorować nimi scatter plot.

## 2. Standaryzacja (autoskalowanie): widget **Preprocess**

Wykład explicit mówi: „**autoskalowanie!**" przed PCA i HCA. Bez tego cechy o dużym zakresie (MW: ~100-1000) zdominują te o małym (FracSP3: 0-1).

- **Preprocess** → dodaj transformację **Normalize Features** → wybierz **Standardize (μ=0, σ²=1)**

> Dla fingerprintów (binarne 0/1) standaryzacja nie jest konieczna — wszystkie cechy mają już ten sam zakres. Możesz Preprocess pominąć dla `combined_morgan.csv` / `combined_maccs.csv` (albo zastosować i sprawdzić różnicę — w prezentacji można o tym wspomnieć).

## 3. PCA: widget **PCA**

Podłącz wyjście `Preprocess → PCA`.

Ustawienia:
- **Components**: zacznij od 10 (zobaczysz scree plot — z reguły PC1+PC2 wyjaśnią 40-70% wariancji deskryptorów ciągłych)
- zaznacz **Show only first n principal components** = 2 lub 3 (do scatter plotu)

**Screenshot 1 do prezentacji** → scree plot pokazujący % wyjaśnianej wariancji per komponent. Z tego liczysz „ile komponentów wystarcza" (regułą kciuka: gdzie się załamuje — kryterium łokcia).

## 4. Wizualizacja: widget **Scatter Plot**

Połącz PCA → Scatter Plot.

Konfiguracja (powtórz w kilku wariantach!):

| Wariant | X | Y | Color | Shape | Label |
|---|---|---|---|---|---|
| A — „grupy" | PC1 | PC2 | `group` | `subgroup` | `name` |
| B — „logP" | PC1 | PC2 | `logP` (continuous) | `group` | `name` |
| C — „MW" | PC1 | PC2 | `MW` (continuous) | `group` | `name` |
| D — „TPSA" | PC1 | PC2 | `TPSA` (continuous) | `group` | `name` |
| E — „głębiej" | PC1 | PC3 | `group` | `subgroup` | `name` |

> **Lifehack**: w Scatter Plot zaznacz opcję „Show class density" — dostajesz miłe konturowe obwiednie wokół klas. Wygląda znakomicie na slajdach.

**Screenshot 2-6 do prezentacji** → po jednym scatter na wariant. Wariant A pokazuje „czy 3 grupy się rozdzielają"; B-D pokazują „o jakim *gradiencie* mówi PC1/PC2".

## 5. HCA — dendrogram: widgety **Distances** + **Hierarchical Clustering**

To jest dokładnie ten workflow, który był na slajdach.

**Distances** ← podłącz z `Preprocess` (nie z PCA!):
- **Distance Metric**: zacznij od *Euclidean* (klasyka z wykładu). Później sprawdź też *Manhattan* — w deskryptorach ciągłych daje czasem ciekawsze podziały.
- **Distance between**: Rows

**Hierarchical Clustering** ← podłącz z Distances:
- **Linkage**: zacznij od **Ward** (zwykle daje najczystsze klastry); później sprawdź **Complete** (najdalszy sąsiad) i **Single** (najbliższy sąsiad — często „daisy-chain" effect).
- **Annotation**: wybierz `name` — będziesz widzieć etykiety na końcach gałęzi
- **Color labels by**: `group`

Wykład explicit zaleca **kryterium Sneatha**: tnij dendrogram na 2/3 wysokości (D). W widget'cie przeciągnij linię cięcia (pojawia się przerywaną linią). Możesz też ustawić "Top N clusters" = 3, 6, 10 i porównać.

**Screenshot 7-9 do prezentacji** → dendrogram dla każdej reprezentacji (deskryptory / Morgan / ChemBERTa). Podpisz każdą gałąź, na której pojawia się ciekawy mostek.

> Dla **Morgan FP** użyj distance metric = *Jaccard* lub *Cosine* zamiast Euclidean. To standard dla binarnych fingerprintów (Tanimoto = Jaccard). Orange nie ma Tanimoto eksplicytnie, ale Jaccard jest matematycznie tym samym dla binarnych danych.

## 6. (Opcjonalnie) k-Means: widget **k-Means**

Workflow: `Preprocess → k-Means → Scatter Plot` (lub `Data Table`).

- **Number of clusters**: ustaw `Fixed: 3` (zgodnie z grupą), potem 6 (podgrupy chemiczne), potem 10
- **Silhouette score** — Orange pokazuje to obok. Im wyższy, tym wyraźniejsze klastry. Zrób tabelkę silhouette dla k=2..15 i wskaż „łokieć".

W Scatter Plot ustaw **Color = Cluster** (z k-Means) i **Shape = group** (z File). Jeśli kształty i kolory pokrywają się dla większości punktów → klastry odpowiadają grupom.

**Screenshot 10 do prezentacji** → silhouette vs k. To jest analiza wrażliwości na liczbę klastrów (punkt 7 briefu).

## 7. (Opcjonalnie) t-SNE: widget **t-SNE**

Pojawia się w niektórych wersjach Orange (Bioinformatics add-on). Jeśli go masz:
- **Perplexity**: testuj 5, 15, 30. Dla 90 punktów *perplexity=15-20* jest sensowne.
- **PCA preprocessing**: zostaw włączone (przyspiesza).

t-SNE jest szczególnie wartościowy dla **fingerprintów** — PCA na binarnych danych słabo działa, t-SNE radzi sobie lepiej.

## 8. Wizualizacja molekuł

Orange ma **Chem add-on**. Jeśli nie chce Ci się instalować, masz folder `mol_images/` — wkleisz obrazki ręcznie na slajdy obok wykresów.

Alternatywnie: workflow `File → Scatter Plot` + ustaw `Show tooltip with image` (Orange wczytuje obrazki ze ścieżki, jeśli kolumna z nazwą pliku jest w danych).

## 9. Pełny canvas Orange (mapa workflow)

```
                                    ┌─→ Scatter Plot (PC1 vs PC2, color=group)
                                    │
File ──→ Preprocess ──→ PCA ────────┼─→ Scatter Plot (PC1 vs PC2, color=logP)
(combined_      (Standardize)       │
 continuous)                        └─→ Scatter Plot (PC1 vs PC3, color=group)
       │
       └─→ Distances (Euclidean) ──→ Hierarchical Clustering (Ward) ─→ Dendrogram
       │
       └─→ k-Means (k=3,6,10) ─────→ Scatter Plot (cluster vs group)
```

Powtórz cały graf dla `combined_morgan.csv` (z Distance = Jaccard) i `combined_chemberta.csv` (Euclidean). **Trzy równoległe pipeline'y = serce porównania reprezentacji.**

## 10. Co konkretnie pokazać na slajdach prezentacji

Slajd „Wyniki — deskryptory ciągłe":
- Scree plot PCA (% wariancji)
- Scatter PC1×PC2 pokolorowany wg `group`
- Scatter PC1×PC2 pokolorowany wg `logP` (continuous) — pokaże gradient hydrofobowości
- Dendrogram HCA (Ward, Euclidean) z linią cięcia w 2/3D

Slajd „Wyniki — Morgan fingerprints":
- Scatter PC1×PC2 (uwaga: % wariancji będzie niski — to ważna obserwacja!)
- Dendrogram HCA (Ward, **Jaccard**)
- t-SNE jeśli masz add-on

Slajd „Wyniki — ChemBERTa embeddings":
- Scatter PC1×PC2
- Dendrogram HCA (Ward, Euclidean lub Cosine)

Slajd „Porównanie 3 reprezentacji" — TRZY scatter ploty obok siebie, ten sam kolor=group. Pokaż, że klastry mają **różne kształty** w zależności od reprezentacji. To jest ten "main result" projektu.

Slajd „Mostki nieoczywiste":
- Zoom na region, gdzie kwas moczowy zbliża się do kwasu askorbinowego (deskryptory ciągłe)
- Zoom na region paracetamol/kofeina
- Pokaż obrazki tych molekuł obok

Slajd „Wrażliwość":
- 3 dendrogramy z różnym linkagiem (Single/Complete/Ward) → ten sam zbiór, inne podziały
- Silhouette vs k z k-Means
- Dla porównania: PCA bez standaryzacji vs z (kolosalna różnica)

## 11. Pułapki, na które uważaj

1. **Zapomniana standaryzacja** → PCA zdominuje MW, scree plot pokaże 95% wariancji w PC1, klastry zdegenerowane. Diagnoza: jeśli loadings PC1 to prawie sam MW/HeavyAtoms — wróć do Preprocess.

2. **Niepokolorowane podgrupy** → na pierwszy rzut oka 3 chmurki, OK. Pokoloruj wg `subgroup` (kontrol-klik) — zobaczysz że np. salicylates i NLPZ_arylpropionic są w innych częściach „chmury Painkillers". To jest *naprawdę* ciekawy obraz.

3. **Fingerprinty + Euclidean distance** → działa, ale to nie najlepszy wybór. Jaccard/Tanimoto to standard branżowy.

4. **Etykiety w Scatter Plot** → przy 90 punktach napisy się nakładają. Pokaż tylko `name` dla *Selected* punktów (kliknięcie/zaznaczenie regionu) — albo użyj Tooltip.

5. **Outliery** → β-karoten i MK-7 (oba ogromne, lipofilowe) prawdopodobnie wylądują daleko na PCA i będą „rozciągać" oś. Pokaż jeden wykres z nimi, drugi po wyrzuceniu — to jest *świetny* przykład analizy wrażliwości do punktu 7 briefu.

6. **„Group" jako feature zamiast jako target/meta** → wtedy PCA „widzi" etykietę i tworzą się klastry przez przypadek. *Group nigdy nie wchodzi do PCA jako feature.*
