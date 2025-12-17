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

