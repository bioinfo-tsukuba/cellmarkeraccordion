# The Cell Marker Accordion 

Marker gene preparation
1. get `TheCellMarkerAccordion_database_v1.0.0.xlsx` from https://github.com/TebaldiLab/shiny_cellmarkeraccordion/tree/main/data 

2. save Human_healthy sheet and Mouse_healthy sheet as `$PATH/cellmarkeraccordion/v.1.0.0/original_human_healthy.csv` and `original_mouse_heaelthy.csv`

3. run `convert_original_to_light.py` in terminal under cell-io-venv 

```sh
conda activate cell-io-venv
cd $PATH/cellmarkeraccordion
python convert_original_to_light.py 

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