---
output:
  pdf_document: default
  html_document: default
---
# Dominant migration strategies within the Gulf of St. Lawrence capelin stock inferred from otolith chemistry transects

**Romaric Jac, Olivier Le Pape, Elisabeth van Beveren, Mathieu Boudreau, Lola Coussau, Pascal Sirois, Dominique Robert, Pablo Brosset**

*Canadian Journal of Fisheries and Aquatic Sciences* — Small Pelagic Fish Symposium (SPF-2026) special issue

> DOI script: https://doi.org/10.5281/zenodo.22112402
> DOI article: *to be added upon publication*

---

## Overview

This repository contains the R script used to produce all analyses and figures presented in the paper. The workflow starts from cleaned otolith chemistry transects and proceeds through:

1. QDA-based region assignment along each transect
2. Multi-resolution binning and sensitivity analyses
3. MCA + k-means clustering of individual migration histories
4. Characterisation of migration strategies and comparison with trawl-survey fork-length data

---

## Repository structure

```
.
├── otolith_transect_analysis.R   # Full analysis script (single file)
├── README.md                     # This file
└── LICENSE                       # MIT licence
```

Output figures and intermediate objects are written to the working directory at runtime. One subfolder is created automatically:

- `plots_cluster_majority_bis/` — bubble plots, one PNG per cluster

---

## Data availability

Otolith chemistry data are **not publicly archived** due to ongoing research programmes. All input files are available upon reasonable request to the corresponding author:

**Romaric Jac** — romaric.jac@uqar.ca

### Input files

| File | Description |
|------|-------------|
| `transect_capelin_clean.csv` | Cleaned LA-ICP-MS transect data. One row per measurement point. |
| `otolith_margin_clean_region.csv` | Otolith edge measurements used to train the QDA model. |
| `FL_nGSL.csv` | Trawl fork-length data, northern Gulf of St. Lawrence. |
| `FL_sGSL.csv` | Trawl fork-length data, southern Gulf of St. Lawrence. |

### Column descriptions

#### `transect_capelin_clean.csv`

| Column | Type | Description |
|--------|------|-------------|
| `Individual` | character | Unique otolith identifier |
| `Distance_to_core` | numeric | Distance from otolith core (µm) |
| `Length_at_dist` | numeric | Reconstructed fork length at that distance (cm) |
| `Length` | numeric | Capture fork length (cm) |
| `Sex` | character | `F`, `M`, or `undetermined` |
| `Survey` | character | `Commercial` or `MPO` |
| `Region_caught` | character | Capture region (`estuary`, `strait`, `south`, `north`) |
| `Region_pred` | character | Region predicted from a prior classification step |
| `Li7_Ca` | numeric | ⁷Li/Ca elemental ratio (µmol mol⁻¹) |
| `B11_Ca` | numeric | ¹¹B/Ca elemental ratio (µmol mol⁻¹) |
| `Mg25_Ca` | numeric | ²⁵Mg/Ca elemental ratio (µmol mol⁻¹) |
| `K39_Ca` | numeric | ³⁹K/Ca elemental ratio (µmol mol⁻¹) |
| `Zn64_Ca` | numeric | ⁶⁴Zn/Ca elemental ratio (µmol mol⁻¹) |
| `Sr88_Ca` | numeric | ⁸⁸Sr/Ca elemental ratio (µmol mol⁻¹) |
| `Ba138_Ca` | numeric | ¹³⁸Ba/Ca elemental ratio (µmol mol⁻¹) |

*Columns 18 and 22–28 in the raw file correspond to the elemental ratio columns above and are averaged during transect binning; all other columns take the first value of each bin.*

#### `otolith_margin_clean_region.csv`

Same elemental ratio columns as above, plus:

| Column | Type | Description |
|--------|------|-------------|
| `Region` | character | Known capture region used as QDA training label |

#### `FL_nGSL.csv`

| Column | Type | Description |
|--------|------|-------------|
| `year` | integer | Survey year |
| `latitude` | numeric | Latitude in DDMM.mm sexagesimal format |
| `longitude` | numeric | Longitude in DDMM.mm sexagesimal format |
| `length` | numeric | Fork length (mm) |
| `number.caught` | numeric | Abundance weight for this length class |

#### `FL_sGSL.csv`

| Column | Type | Description |
|--------|------|-------------|
| `year` | integer | Survey year |
| `length` | numeric | Fork length (cm) |
| `number.caught` | numeric | Abundance weight for this length class |

---

## Figures produced

| File | Section | Description |
|------|---------|-------------|
| `correlation.png` | §4 | Pairwise correlations among elemental ratios |
| `State_distribution_plot.png` | §6 | Region proportions along the fork-length axis |
| `Transition_diagram.png` | §7 | Directed graph of permanent habitat transitions |
| `Barplot_nombre_transitions.png` | §7 | Distribution of transition counts per individual |
| `Total_transitions_vs_threshold.png` | §8 | Sensitivity: total transitions vs detection threshold |
| `Proportion_vs_threshold.png` | §8 | Sensitivity: mobility-class proportions vs threshold |
| `undetermined.png` | §9 | Histogram of unidentified-region proportions |
| `variances.png` | §10 | MCA scree plot |
| `Clusters.png` | §11 | k-means clusters in MCA space |
| `Barplot_clusters.png` | §11 | Capture region × survey per cluster |
| `Barplot_clusters_sex.png` | §11 | Sex composition per cluster |
| `plots_cluster_majority_bis/cluster_*_bubbles.png` | §12 | Synthetic profiles per cluster |
| `shift_SBI_sud.png` | §13 | SBI vs South proportions along length (cluster 3) |
| `density_NGSL_SGSL.png` | §13 | Weighted fork-length distributions, n-GSL vs s-GSL |

---

## How to run

1. Clone or download this repository.
2. Place the four input CSV files in the **same directory** as the script (or set your R working directory accordingly).
3. Open R (≥ 4.3) or RStudio and run:

```r
source("otolith_transect_analysis.R")
```

All figures are saved to the working directory automatically. Expected runtime is a few minutes depending on dataset size.

---

## Dependencies

Install all required packages before running:

```r
install.packages(c(
  "MASS", "dplyr", "tidyr", "ggplot2", "purrr", "readr",
  "cluster", "factoextra", "FactoMineR",
  "fpc", "ggsci", "colorspace", "ggnewscale", "ggforce",
  "nnet", "patchwork", "scales"
))
```

The script was developed and tested under the following environment. Full details are available in [`session_info.txt`](session_info.txt).

| Component | Version |
|-----------|---------|
| R | 4.4.2 (2024-10-31 ucrt) |
| Platform | x86\_64-w64-mingw32/x64 (Windows 11, build 26100) |
| **Attached packages** | |
| MASS | 7.3-61 |
| ggplot2 | 4.0.2 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |
| purrr | 1.2.1 |
| readr | 2.1.5 |
| cluster | 2.1.8.2 |
| factoextra | 1.0.7 |
| FactoMineR | 2.13 |
| fpc | 2.2-14 |
| ggsci | 3.2.0 |
| colorspace | 2.1-1 |
| ggnewscale | 0.5.2 |
| ggforce | 0.4.2 |
| nnet | 7.3-19 |
| patchwork | 1.3.0 |
| scales | 1.4.0 |
| Matrix | 1.7-1 |
| lme4 | 1.1-35.5 |

---

## Licence

This code is released under the MIT Licence — see [`LICENSE`](LICENSE) for details.
