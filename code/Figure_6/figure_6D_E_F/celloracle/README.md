# CellOracle analysis for HIF1A knockout

Reproducible CellOracle workflow associated with:

> **HIF-1α integrates metabolic and immunoregulatory programs in RORγt⁺ regulatory T cells during intestinal inflammation**
>
> Cipelli et al.

## What generates what

**[Open the complete pipeline guide as a PDF table](output/pdf/CellOracle_HIF1A_pipeline_guide.pdf)**

| Order | Script | Reads | Generates |
|---:|---|---|---|
| 1 | `run_hif1a_knockout.py` | Seurat export: counts, genes, cells, metadata and UMAP | Fitted cluster-specific GRNs, simulated HIF1A-KO Oracle, Links, delta tables and primary UMAP/vector figures |
| 2 | `generate_FA_paperstyle_HIF1A_KO.py` | Simulated Oracle and cell-level delta metrics | ForceAtlas embeddings, cell/grid vectors, L2 magnitude tables and publication-ready FA figures |
| 3 | `generate_PS_markov_HIF1A_KO.py` | Simulated Oracle with UMAP, pseudotime and group metadata | Pseudotime gradient, perturbation score, Markov trajectories, endpoint-density tables and PS/Markov figures |
| 4 | `regenerate_propagation_HIF1A_KO.py` | Existing simulated Oracle | Final L1 `Mean delta X length` propagation table and PNG/PDF figure |

Scripts 2–4 are independent downstream branches. Each uses the Oracle produced by script 1.

```text
Seurat export
    |
    v
run_hif1a_knockout.py
    |
    +--> ForceAtlas figures
    +--> perturbation score + Markov
    +--> propagation-impact figure
```

For installation, complete commands, output filenames, scientific interpretation and the distinction between L1 and L2 metrics, see [README_hif1a_knockout.md](README_hif1a_knockout.md).
