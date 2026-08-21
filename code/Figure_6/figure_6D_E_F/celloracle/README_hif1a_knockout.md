# CellOracle HIF1A Knockout

Pipeline para simular nocaute in silico de `HIF1A` no objeto FOXP3+ Treg
`data/seu_obj/SC_obj_sub_filtered.rds`.

## 1. Exportar Seurat para CellOracle

Execute a partir de `prof_niels`:

```bash
Rscript code/knowkout.r
```

Saida:

```text
celloracle_export_hif1a/
  counts_genes_x_cells.mtx
  genes.tsv
  cells.tsv
  metadata.csv
  umap.csv
  export_summary.csv
```

O export usa `RNA/counts`, conserva `seurat_clusters`, `group`, UMAP e, se existir,
`data/seu_obj/monocle2_cds_ordered.rds`, adiciona `Pseudotime` e `State`.
Por padrao exporta 3.000 genes, mantendo `HIF1A`.

Parametros por variavel de ambiente:

```bash
CELLO_TARGET_GENE=HIF1A \
CELLO_N_FEATURES=3000 \
CELLO_SEURAT_RDS=data/seu_obj/SC_obj_sub_filtered.rds \
CELLO_EXPORT_DIR=celloracle_export_hif1a \
Rscript code/knowkout.r
```

## 2. Preparar Python

CellOracle nao estava instalado no ambiente Python atual. O Python local era 3.13,
mas `celloracle==0.20.0` prende `pandas<=1.5.3` e `matplotlib<3.7`, entao use
Python 3.10 ou 3.11.

Opcao recomendada com conda/mamba:

```bash
mamba env create -f code/Figure_6/figure_6D_E_F/celloracle/environment.yml
conda activate celloracle-hif1a
bash code/Figure_6/figure_6D_E_F/celloracle/install_celloracle_extras.sh celloracle-hif1a
```

Opcao com `venv`, se `python3.10` estiver instalado:

```bash
python3.10 -m venv .venv-celloracle
source .venv-celloracle/bin/activate
python -m pip install -r code/Figure_6/figure_6D_E_F/celloracle/requirements-celloracle.txt
```

No macOS arm64 deste projeto, `velocyto` precisou ser compilado com
`code/Figure_6/figure_6D_E_F/celloracle/clang_openmp_wrapper.sh`, pois o Apple clang nao aceita
`-fopenmp` diretamente.

## 3. Rodar o nocaute in silico

```bash
python code/Figure_6/figure_6D_E_F/celloracle/run_hif1a_knockout.py \
  --export-dir celloracle_export_hif1a \
  --out-dir results/celloracle_hif1a_ko \
  --target-gene HIF1A \
  --species human \
  --genome hg38 \
  --cluster-column seurat_clusters \
  --group-column group
```

Rodada concluida neste workspace, em modo mais leve por causa das restricoes de
multiprocessing do sandbox:

```bash
CELLO_N_FEATURES=1500 CELLO_EXPORT_DIR=celloracle_export_hif1a_1500 \
  Rscript code/knowkout.r

conda run -n celloracle-hif1a python code/Figure_6/figure_6D_E_F/celloracle/run_hif1a_knockout.py \
  --export-dir celloracle_export_hif1a_1500 \
  --out-dir results/celloracle_hif1a_ko_1500 \
  --target-gene HIF1A \
  --species human \
  --genome hg38 \
  --cluster-column seurat_clusters \
  --group-column group \
  --n-jobs 1 \
  --bagging-number 10 \
  --link-threshold-number 1000
```

O script usa o pipeline principal do CellOracle:

1. cria `Oracle` a partir de `AnnData`;
2. importa promoter base-GRN humana;
3. roda PCA e KNN imputation;
4. infere `Links` por subcluster;
5. filtra links;
6. ajusta GRN cluster-especifica;
7. simula `perturb_condition={"HIF1A": 0.0}`;
8. calcula probabilidades de transicao e vetores no UMAP.

## Objetos e tabelas principais

```text
results/celloracle_hif1a_ko/
  objects/
    input_for_celloracle.h5ad
    oracle_imported.celloracle.oracle
    oracle_preprocessed.celloracle.oracle
    links_filtered.celloracle.links
    oracle_fitted_grn.celloracle.oracle
    oracle_HIF1A_KO.celloracle.oracle
  tables/
    base_grn_used.csv
    links_raw.csv
    links_filtered.csv
    network_scores.csv
    gene_delta_summary.csv
    cell_delta_metrics.csv
    cluster_delta_summary.csv
  figures/
    umap_cluster_and_target.png
    hif1a_ko_quiver_cells.png
    hif1a_ko_grid_arrows.png
    n_propagation_impact.png
  logs/
    run_info.json
```

Na rodada concluida, os artefatos equivalentes estao em:

```text
results/celloracle_hif1a_ko_1500/
celloracle_export_hif1a_1500/
```

Se voce tiver uma base GRN propria derivada de scATAC/motif analysis, use:

```bash
python code/Figure_6/figure_6D_E_F/celloracle/run_hif1a_knockout.py \
  --base-grn path/to/TF_info_matrix.parquet
```

## Organizacao dos scripts Python

Os arquivos executaveis mantem apenas argumentos e etapas de pipeline. As funcoes
reutilizaveis ficam em `functions/`:

- `io_utils.py`: ambiente, diretorios e leitura das tabelas de entrada;
- `celloracle_utils.py`: carga/validacao de objetos, GRN e exportacao de tabelas;
- `embedding_utils.py`: ForceAtlas, embeddings, vetores, projecoes e grid flow;
- `perturbation_utils.py`: propagation magnitude, delta L2, PS e Markov;
- `plotting_utils.py`: plots, exportacao PNG/PDF e a paleta `CLUSTER_COLORS`.

O grafico publication-ready de propagation magnitude e gerado por
`plot_propagation_magnitude()` em `functions/plotting_utils.py`, sem normalizar as
curvas e usando as cores fixas dos clusters 0--6.

Para regenerar apenas a figura original `Mean delta_X length` a partir do Oracle
ja calculado, sem repetir o pipeline completo:

```bash
conda run -n celloracle-hif1a python \
  code/Figure_6/figure_6D_E_F/celloracle/regenerate_propagation_HIF1A_KO.py
```
