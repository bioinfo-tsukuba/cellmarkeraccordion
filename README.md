# The Cell Marker Accordion 

Marker gene preparation
1. get `TheCellMarkerAccordion_database_v1.0.0.xlsx` from https://github.com/TebaldiLab/shiny_cellmarkeraccordion/tree/main/data 

2. save Human_healthy sheet and Mouse_healthy sheet as `$PATH/cellmarkeraccordion/v.1.0.0/original_human_healthy.csv` and `original_mouse_heaelthy.csv`

3. to add ensemble id, it uses biomart. to change release version, modify `biomart_config.json`

```
{
  "marker_list_set": [
    {
      "id": 1,
      "list": [
        {
          "version": "1.0.0",
          "biomart_configs": [
            {
              "role": "source",
              "file": "original_human_healthy.csv",
              "outfile": "light_human_healthy.csv",
              "species": "Human",
              "tax_id": 9606,
              "ensembl_dataset": "hsapiens_gene_ensembl",
              "symbol_attribute": "hgnc_symbol",
              "backend": "ensembl_biomart",
              "release": 110
            },
            {
              "role": "source",
              "file": "original_mouse_healthy.csv",
              "outfile": "light_mouse_healthy.csv",
              "species": "Mouse",
              "tax_id": 10090,
              "ensembl_dataset": "mmusculus_gene_ensembl",
              "symbol_attribute": "mgi_symbol",
              "backend": "ensembl_biomart",
              "release": 110
            },
            {
              "role": "derived",
              "outfile": "light_rat_healthy.csv",
              "species": "Rat",
              "tax_id": 10116,
              "ensembl_dataset": "rnorvegicus_gene_ensembl",
              "symbol_attribute": "rgd_symbol",
              "source_species": "Mouse",
              "backend": "ensembl_biomart",
              "release": 110
            }
          ]
        }
      ]
    }
  ]
}
```

4. run `convert_original_to_light.py` in terminal under cell-io-venv 

```sh
conda activate cell-io-venv # set up via /data01/iharuto/cell-io-mappingctl/0_VENV_SETUP.sh or https://github.com/bioinfo-tsukuba/cell-io-mappingctl/blob/main/0_VENV_SETUP.sh
cd $PATH/cellmarkeraccordion
python convert_original_to_light.py # for more config id, python convert_original_to_light.py --config-id 2

2025-12-15 16:20:55 - Python.conversion_pipeline - INFO - Converting markers for config ID 1
2025-12-15 16:20:55 - Python.conversion_pipeline - INFO - Processing Human...
2025-12-15 16:21:01 - Python.conversion_pipeline - INFO - Saved Human: 171052 records
2025-12-15 16:21:01 - Python.conversion_pipeline - INFO - Processing Mouse...
2025-12-15 16:21:07 - Python.conversion_pipeline - INFO - Saved Mouse: 49419 records
2025-12-15 16:21:07 - Python.conversion_pipeline - INFO - Processing Rat (from Mouse)...
2025-12-15 16:21:10 - Python.conversion_pipeline - INFO - Saved Rat: 56222 records
2025-12-15 16:21:10 - Python.conversion_pipeline - INFO - Conversion completed!
2025-12-15 16:21:10 - __main__ - INFO - Conversion completed!
```

## Cell Type Scoring Methodology

CellMarkerAccordion uses a weighted scoring system to assign cell types to clusters based on marker gene expression. The scoring combines expression levels with evidence quality and marker specificity.

### Score Components

#### **ECs_global (Evidence Consistency Score)**
- **Definition**: Counts how many different studies/resources report the same marker for the same cell type
- **Calculation**: `ECs_global = number_of_independent_sources`
- **Range**: 1 to 15+ (higher = better)
- **Interpretation**: 
  - ECs = 1: Only one study reports this marker for this cell type
  - ECs = 7: Seven different studies independently report this marker
  - ECs = 15: Fifteen studies report this marker (very high confidence)

#### **SPs_global (Specificity Score)**
- **Definition**: Measures how specific a marker is to particular cell types
- **Calculation**: `SPs_global = 1 / number_of_cell_types_using_this_marker`
- **Range**: 0.01 to 1.00 (higher = better)
- **Interpretation**:
  - SPs = 1.00: Perfectly specific (only 1 cell type uses this marker)
  - SPs = 0.50: Moderate specificity (2 cell types use this marker)
  - SPs = 0.06: Low specificity (~17 cell types use this marker)

### Cell Type Score Calculation

For each cluster and each candidate cell type, the algorithm calculates a weighted score:

```python
# For each cell type
total_score = 0.0
total_weight = 0.0

for marker in positive_markers_for_cell_type:
    expr_level = cluster_mean_expression[marker]  # Mean expression in cluster
    ecs_score = marker_ecs_global                 # Evidence score
    sps_score = marker_sps_global                 # Specificity score
    
    weight = ecs_score * sps_score
    weighted_expr = expr_level * weight
    
    total_score += weighted_expr
    total_weight += weight

cell_type_score = total_score / total_weight  # Weighted average
```

### Example Calculation

For a "T cell" annotation with markers CD3E, CD8A, CD4:

| Marker | Expression | ECs | SPs | Weight | Weighted_Expr |
|--------|------------|-----|-----|--------|---------------|
| CD3E   | 5.2        | 12  | 0.8 | 9.6    | 49.92        |
| CD8A   | 3.1        | 8   | 0.9 | 7.2    | 22.32        |
| CD4    | 1.8        | 10  | 0.7 | 7.0    | 12.60        |

```
total_score = 49.92 + 22.32 + 12.60 = 84.84
total_weight = 9.6 + 7.2 + 7.0 = 23.8
T_cell_score = 84.84 / 23.8 = 3.56
```

### Final Assignment

The cluster is assigned to the cell type with the highest weighted score. This approach:
- **Balances evidence and specificity**: High-evidence but non-specific markers don't dominate
- **Accounts for expression levels**: Actually expressed markers contribute more
- **Provides robust annotation**: Multiple independent evidence sources increase confidence

