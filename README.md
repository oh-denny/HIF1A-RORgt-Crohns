# HIF-1α integrates metabolic and immunoregulatory programs in RORγt⁺ regulatory T cells during intestinal inflammation

Code and supporting documentation for the analyses reported in the manuscript:

> **HIF-1α integrates metabolic and immunoregulatory programs in RORγt⁺ regulatory T cells during intestinal inflammation**<br>
> Marcella Cipelli*, Eloísa M. da Silva, Luísa Menezes-Silva, Barbara N. Padovani, Mariana A. Amaral, Laís C. Paredes, Bruno G. Nunes, Victor Y. Yariwake, José Arimatéia O. N. Neto, Natalia N. Bos, Anthony G. da Silveira, João Vinicius H. da Silva, Raquel S. Vieira, Suemy M. Yamada, Luis Felipe S. Moreira, Benedito Matheus dos Santos, Aline Ignacio, Maria Fernanda Forni, Orestes Foresto-Neto, Jefferson Antonio Leite, Marco Aurélio R. Vinolo, Dennyson Leandro M Fonseca, Sandra Marcia Muxel, Matthias Lochner, Vinicius Andrade-Oliveira, and Niels O. S. Camara*.

`*` Corresponding authors.

## Overview

This repository contains the computational workflow used to reanalyze human ileal single-cell RNA-sequencing data from Crohn's disease and healthy samples. The analyses focus on FOXP3⁺ T-cell states, RORγt- and HIF-1α-associated transcriptional programs, pseudotime trajectories, and CellOracle simulations of an *HIF1A* in silico knockout.

The main workflow is:

```text
GSE209832 data
      ↓
Seurat quality control and SCTransform normalization
      ↓
FOXP3⁺ T-cell subclustering and transcriptional analyses
      ├── Figure 1: cell states, HIF1A, RORγt, and pathway enrichment
      └── Figure 6: Monocle trajectory and pathway scores
                     ↓
              CellOracle GRN analysis
                     ↓
          HIF1A knockout simulation (Figures 6D–F)
```

## Repository structure

| Path | Purpose | Main output |
|---|---|---|
| [`code/data_download`](code/data_download) | Downloads and organizes the GSE209832 single-cell dataset | Raw HDF5 matrices and sample metadata |
| [`code/processing`](code/processing) | Seurat quality control, SCTransform normalization, integration, and clustering | Processed Seurat objects |
| [`code/Figure_1`](code/Figure_1) | FOXP3⁺ subclustering, RORγt/HIF1A analyses, pseudobulk differential expression, and enrichment | Figure 1 panels and result tables |
| [`code/Figure_6`](code/Figure_6) | Monocle trajectory, pathway scores, expression heatmap, and perturbation analysis | Figure 6 panels |
| [`code/Figure_6/figure_6D_E_F/celloracle`](code/Figure_6/figure_6D_E_F/celloracle) | CellOracle network construction and *HIF1A* knockout simulation | Figures 6D–F and intermediate CellOracle objects |
| [`code/Supp`](code/Supp) | Supplementary RORγt-high abundance and donor/sample analyses | Supplementary figures |
| [`functions`](functions) | Shared R functions and Monocle compatibility patches | Reusable analysis helpers |
| [`figures`](figures) | Generated figure files | TIFF and related figure outputs |
| [`pseudobulk_clusters`](pseudobulk_clusters) | Pseudobulk differential-expression and enrichment results | Tables and pathway plots |

## What generates each figure

| Figure or analysis | Script | What it generates |
|---|---|---|
| Figure 1A | [`code/Figure_1/figure_1A/figure_1A.r`](code/Figure_1/figure_1A/figure_1A.r) | UMAP of FOXP3⁺ T-cell subclusters |
| Figure 1B | [`code/Figure_1/Figure_1B/figure_1B.r`](code/Figure_1/Figure_1B/figure_1B.r) | Distribution of RORγt-high Tregs by study group |
| Figure 1C | [`code/Figure_1/figure_1C/figure_1C.r`](code/Figure_1/figure_1C/figure_1C.r) | Split violin plots of Treg- and HIF-related markers |
| Figure 1D | [`code/Figure_1/figure_1D/figure_1D.r`](code/Figure_1/figure_1D/figure_1D.r) | Alluvial plot connecting group, FOXP3⁺ subcluster, and HIF1A status |
| Figure 1E | [`code/Figure_1/figure_1E/figure_1E.r`](code/Figure_1/figure_1E/figure_1E.r) | Donut charts showing the distribution of HIF1A⁺ cells |
| Figure 1F | [`code/Figure_1/figure_1F/figure_1F.r`](code/Figure_1/figure_1F/figure_1F.r) | Pseudobulk ssGSEA heatmap for selected clusters |
| Figure 1G | [`code/Figure_1/figure_1G/figure_1G.r`](code/Figure_1/figure_1G/figure_1G.r) | Pseudobulk DESeq2 and Hallmark fgsea analyses, including hypoxia enrichment |
| Figure 6 trajectory input | [`code/Figure_6/trajectory.r`](code/Figure_6/trajectory.r) | Ordered Monocle 2 CellDataSet used by the trajectory panels |
| Figure 6A | [`code/Figure_6/figure_6A/figura_6A.r`](code/Figure_6/figure_6A/figura_6A.r) | Trajectory colored by pseudotime and Monocle state |
| Figure 6B | [`code/Figure_6/figure_6B/figure_6B.r`](code/Figure_6/figure_6B/figure_6B.r) | Trajectory colored by cell cluster and clinical group |
| Figure 6C | [`code/Figure_6/figure_6C/figure_6C.r`](code/Figure_6/figure_6C/figure_6C.r) | Hypoxia and TNFα/NF-κB UCell scores along the trajectory |
| Figures 6D–F | [`code/Figure_6/figure_6D_E_F/celloracle`](code/Figure_6/figure_6D_E_F/celloracle) | CellOracle GRNs, *HIF1A* knockout propagation, vector fields, and perturbation scores |
| Figure 6G | [`code/Figure_6/figure_6G/figure_6G.r`](code/Figure_6/figure_6G/figure_6G.r) | Cluster-level heatmap of Treg, inflammatory, hypoxia, and mitochondrial genes |
| Supplementary analyses | [`code/Supp/Supp1.r`](code/Supp/Supp1.r) | RORγt-high abundance, fold-change, and donor/sample plots |

For a compact, printable map of the CellOracle workflow, see the [CellOracle HIF1A pipeline guide](code/Figure_6/figure_6D_E_F/celloracle/output/pdf/CellOracle_HIF1A_pipeline_guide.pdf). Detailed execution instructions are provided in the [CellOracle README](code/Figure_6/figure_6D_E_F/celloracle/README.md).

## Data source

The human ileal single-cell RNA-sequencing data are publicly available from the NCBI Gene Expression Omnibus under accession [GSE209832](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209832).

The download script uses `GEOquery` to retrieve the supplementary files and construct initial Seurat objects:

```text
code/data_download/download.r
```

Large raw data files, processed objects, and analysis intermediates may not be stored in Git. Their expected locations are referenced by the scripts under `data/` and by the CellOracle documentation.

## Software requirements

The R workflow uses packages including:

- Seurat and sctransform
- GEOquery
- dplyr, tidyr, ggplot2, and ggalluvial
- monocle (Monocle 2)
- UCell
- DESeq2 and fgsea
- GSVA
- ComplexHeatmap and circlize

The Python requirements for the perturbation workflow are documented separately in [`code/Figure_6/figure_6D_E_F/celloracle/environment.yml`](code/Figure_6/figure_6D_E_F/celloracle/environment.yml) and [`requirements-celloracle.txt`](code/Figure_6/figure_6D_E_F/celloracle/requirements-celloracle.txt).

## Reproducibility notes

- Clone the repository and retain its directory structure because the scripts use project-relative paths.
- Most R scripts should be launched from the repository root. The data-download script contains a project-specific working-directory instruction that may need to be adjusted for the local clone.
- Individual figure scripts are analysis-stage scripts. Some expect precomputed objects—such as Seurat, pseudobulk, or Monocle objects—to exist on disk or already be loaded in the R session.
- Generate the ordered Monocle object with `code/Figure_6/trajectory.r` before running the Figure 6 trajectory panels.
- Follow the dedicated CellOracle README for the required Python environment and execution order of Figures 6D–F.
- Exact random seeds, filtering thresholds, cluster selections, and plotting parameters are defined in the corresponding scripts.

## Interpretation of the CellOracle simulation

In the in silico knockout, *HIF1A* expression is set to zero and the predicted expression shift is propagated through cluster-specific gene-regulatory networks. The simulated vector field is compared locally with the pseudotime-gradient vector field. Positive inner-product scores indicate alignment with the observed trajectory, whereas negative scores indicate opposition.

These simulations predict changes in transcriptional state and cell identity; they do not estimate post-perturbation cell numbers.

## Affiliations

1. Department of Immunology, Institute of Biomedical Sciences, University of São Paulo, São Paulo, Brazil.
2. Renal Division, Department of Clinical Medicine, Faculty of Medicine, University of São Paulo, São Paulo, Brazil.
3. Department of Immunology, Institute of Microbiology Paulo de Góes, Center of Health Sciences, Federal University of Rio de Janeiro, Rio de Janeiro, Brazil.
4. Department of Microbiology, Immunology and Infectious Diseases, Snyder Institute for Chronic Diseases, Cumming School of Medicine, University of Calgary, Alberta, Canada.
5. Department of Genetics, Evolution, Microbiology and Immunology, Institute of Biology, University of Campinas, Campinas, Brazil.
6. School of Arts, Sciences and Humanities, University of São Paulo, São Paulo, Brazil.
7. Institute of Medical Microbiology and Hospital Epidemiology, Hannover Medical School, Hannover, Germany.
8. Center for Natural and Human Sciences, Federal University of ABC, São Paulo, Brazil.

## Citation

The manuscript citation and persistent identifier will be added after publication. Until then, please cite the manuscript title and authors shown above when referring to this repository.
