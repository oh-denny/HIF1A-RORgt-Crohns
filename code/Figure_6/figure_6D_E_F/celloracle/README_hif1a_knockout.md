# CellOracle: HIF1A knockout em RORγt⁺ Tregs

Este diretório contém o pipeline CellOracle usado no trabalho:

> **HIF-1α integrates metabolic and immunoregulatory programs in RORγt⁺ regulatory T cells during intestinal inflammation**
>
> Marcella Cipelli, Eloísa M. da Silva, Luísa Menezes-Silva, Barbara N. Padovani, Mariana A. Amaral, Laís C. Paredes, Bruno G. Nunes, Victor Y. Yariwake, José Arimatéia O. N. Neto, Natalia N. Bos, Anthony G. da Silveira, João Vinicius H. da Silva, Raquel S. Vieira, Suemy M. Yamada, Luis Felipe S. Moreira, Benedito Matheus dos Santos, Aline Ignacio, Maria Fernanda Forni, Orestes Foresto-Neto, Jefferson Antonio Leite, Marco Aurélio R. Vinolo, Dennyson Leandro M. da Fonseca, Sandra Marcia Muxel, Matthias Lochner, Vinicius Andrade-Oliveira e Niels O. S. Camara.

## Guia visual em PDF

O fluxo completo, em formato de tabela, está disponível em:

**[Abrir o guia “o que gera o quê”](output/pdf/CellOracle_HIF1A_pipeline_guide.pdf)**

## Visão rápida: o que gera o quê

| Ordem | Script | Entrada principal | O que faz | Principais saídas |
|---:|---|---|---|---|
| 1 | `run_hif1a_knockout.py` | Export Seurat (`counts`, genes, células, metadados e UMAP) | Cria o Oracle, infere e filtra GRNs por cluster, define `HIF1A = 0`, simula o knockout e calcula vetores no UMAP | Oracle simulado, Links, tabelas de delta, UMAP, quiver, grid flow e propagation impact |
| 2 | `generate_FA_paperstyle_HIF1A_KO.py` | Oracle simulado + `cell_delta_metrics.csv` | Recalcula a visualização no ForceAtlas e organiza magnitude, vetores e grid flow | Figuras FA em PNG/PDF, embeddings e vetores em CSV, Oracle FA |
| 3 | `generate_PS_markov_HIF1A_KO.py` | Oracle simulado com UMAP e pseudotempo | Calcula gradiente de pseudotempo, perturbation score e simulação de Markov | Painéis PS/Markov, tabelas por célula/cluster, transições e gradiente CellOracle |
| 4 | `regenerate_propagation_HIF1A_KO.py` | Oracle simulado já existente | Recalcula somente os passos 0–5 do propagation impact, sem repetir GRN/PCA | Figura final `Mean delta X length` em PNG/PDF e tabela correspondente |

Os scripts 2–4 dependem do Oracle produzido pelo script 1 e podem ser executados separadamente depois que o Oracle principal existe.

## Fluxo recomendado

```text
Seurat export
    |
    v
run_hif1a_knockout.py
    |
    +--> generate_FA_paperstyle_HIF1A_KO.py  --> mapas e vetores ForceAtlas
    |
    +--> generate_PS_markov_HIF1A_KO.py      --> perturbation score + Markov
    |
    +--> regenerate_propagation_HIF1A_KO.py  --> propagation impact final
```

## 1. Preparar o ambiente

Use Python 3.10 ou 3.11. O ambiente reproduzível está em `environment.yml`.

```bash
mamba env create -f code/Figure_6/figure_6D_E_F/celloracle/environment.yml
conda activate celloracle-hif1a
bash code/Figure_6/figure_6D_E_F/celloracle/install_celloracle_extras.sh celloracle-hif1a
```

Alternativa com `venv`:

```bash
python3.10 -m venv .venv-celloracle
source .venv-celloracle/bin/activate
python -m pip install -r code/Figure_6/figure_6D_E_F/celloracle/requirements-celloracle.txt
```

No macOS ARM, `clang_openmp_wrapper.sh` contém o suporte usado para compilar dependências que requerem OpenMP.

## 2. Exportar o objeto Seurat

O pipeline Python espera uma pasta contendo:

```text
celloracle_export_hif1a_1500/
├── counts_genes_x_cells.mtx
├── genes.tsv
├── cells.tsv
├── metadata.csv
└── umap.csv
```

Os metadados devem incluir `seurat_clusters` e `group`. Para perturbation score, também devem incluir a coluna de pseudotempo, usada como `Pseudotime` por padrão.

## 3. Executar o knockout principal

```bash
conda run -n celloracle-hif1a python \
  code/Figure_6/figure_6D_E_F/celloracle/run_hif1a_knockout.py \
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

Saídas principais:

```text
results/celloracle_hif1a_ko_1500/
├── objects/
│   ├── input_for_celloracle.h5ad
│   ├── oracle_imported.celloracle.oracle
│   ├── oracle_preprocessed.celloracle.oracle
│   ├── links_filtered.celloracle.links
│   ├── oracle_fitted_grn.celloracle.oracle
│   └── oracle_HIF1A_KO.celloracle.oracle
├── tables/
│   ├── base_grn_used.csv
│   ├── links_raw.csv
│   ├── links_filtered.csv
│   ├── network_scores.csv
│   ├── gene_delta_summary.csv
│   ├── cell_delta_metrics.csv
│   └── cluster_delta_summary.csv
├── figures/
│   ├── umap_cluster_and_target.png
│   ├── hif1a_ko_quiver_cells.png
│   ├── hif1a_ko_grid_arrows.png
│   └── n_propagation_impact.png/.pdf
└── logs/run_info.json
```

## 4. Gerar as figuras ForceAtlas

```bash
conda run -n celloracle-hif1a python \
  code/Figure_6/figure_6D_E_F/celloracle/generate_FA_paperstyle_HIF1A_KO.py
```

Esse script gera mapas de clusters, magnitude de `delta_l2`, vetores por célula, fluxo no grid, distribuições por cluster/grupo e a curva de propagation baseada em L2.

## 5. Gerar perturbation score e Markov

```bash
conda run -n celloracle-hif1a python \
  code/Figure_6/figure_6D_E_F/celloracle/generate_PS_markov_HIF1A_KO.py
```

O perturbation score é o produto interno local entre o fluxo simulado do knockout e o gradiente de pseudotempo. Valores positivos indicam alinhamento com a progressão observada; valores negativos indicam oposição.

## 6. Regenerar somente o gráfico de propagation

```bash
conda run -n celloracle-hif1a python \
  code/Figure_6/figure_6D_E_F/celloracle/regenerate_propagation_HIF1A_KO.py
```

Saídas:

- `results/celloracle_hif1a_ko_1500/figures/HIF1A_KO_propagation_magnitude_by_cluster.png`
- `results/celloracle_hif1a_ko_1500/figures/HIF1A_KO_propagation_magnitude_by_cluster.pdf`
- `results/celloracle_hif1a_ko_1500/tables/HIF1A_KO_propagation_mean_delta_x_length_by_cluster.csv`

## Atenção: L1 e L2 são métricas diferentes

| Nome no pipeline | Norma | Onde aparece |
|---|---|---|
| `Mean delta X length` | L1, padrão do método CellOracle usado para propagation impact | `regenerate_propagation_HIF1A_KO.py` e `n_propagation_impact` |
| `mean_l2_norm` / `delta_l2` | L2, distância Euclidiana no espaço de expressão | Tabelas e figuras produzidas por `generate_FA_paperstyle_HIF1A_KO.py` |

Não descreva essas duas medidas como equivalentes no manuscrito.

## Interpretação científica

- O knockout é simulado definindo a expressão de `HIF1A` como zero.
- A alteração é propagada pelas GRNs específicas de cluster.
- Os vetores no embedding representam direções previstas de mudança de estado celular.
- O perturbation score mede alinhamento ou oposição à progressão de pseudotempo.
- As análises não estimam diretamente número de células, proliferação, sobrevivência ou abundância pós-perturbação.

## Módulos reutilizáveis

| Módulo | Responsabilidade |
|---|---|
| `functions/io_utils.py` | Diretórios, cache e organização de outputs |
| `functions/celloracle_utils.py` | Leitura e validação de Oracle, GRN e exportação de tabelas |
| `functions/embedding_utils.py` | ForceAtlas, embeddings e vetores de deslocamento/grid |
| `functions/perturbation_utils.py` | Propagation, magnitude, perturbation score e Markov |
| `functions/plotting_utils.py` | Figuras publication-ready, exportação PNG/PDF e `CLUSTER_COLORS` |

## Autoria e afiliações

<details>
<summary>Instituições participantes</summary>

1. Department of Immunology, Institute of Biomedical Sciences, University of São Paulo, São Paulo, Brazil.
2. Renal Division, Department of Clinical Medicine, Faculty of Medicine, University of São Paulo, São Paulo, Brazil.
3. Department of Immunology, Institute of Microbiology Paulo de Góes, Center of Health Sciences, Federal University of Rio de Janeiro, Rio de Janeiro, Brazil.
4. Department of Microbiology, Immunology and Infectious Diseases, Snyder Institute for Chronic Diseases, Cumming School of Medicine, University of Calgary, Alberta, Canada.
5. Department of Genetics, Evolution, Microbiology and Immunology, Institute of Biology, University of Campinas, Campinas, Brazil.
6. School of Arts, Sciences and Humanities, University of São Paulo, São Paulo, Brazil.
7. Institute of Medical Microbiology and Hospital Epidemiology, Hannover Medical School, Hannover, Germany.
8. Center for Natural and Human Sciences, Federal University of ABC, São Paulo, Brazil.

</details>
