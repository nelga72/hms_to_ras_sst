Contributors: Marissela Gomez (marissela.gomez@stantec.com), Nada Elgamal (nada.elgamal@stantec.com)
Last updated: 4/1/2026

Folder Contents: 
* hdf_flow_comparison.py
* hdf_flow_comparison_functions.py
* hdf_flow_comparison_nb.ipynb
* ReadMe.txt

To run: 
* Be sure to activate the environment with necessary packages either in terminal or your code editor (regularhome.yml)

In terminal: 
* Change directory to where the code lives: cd [filepath]. 
** If code is on a different drive (e.g. D:/ drive) change drives before navigating to directory using (cd d:)
* Once in code directory, type "python hdf_flow_comparison.py" and hit Enter.  
* Once prompted, enter the folder path with the model outputs to compare, this folder path should end with the HUC10 number (e.g. D:\Folder1\Folder2\models\1008001203) and hit Enter. 

In Jupyter Notebook:
* Open hdf_flow_comparison_nb.ipynb
* Replace the "inputs_filepath" variable with the filepath of the HUC10 folder.

Description:
This script reviews the reference line flow in the p##.hdf file and input flow from the dss hydrograph in the Hydrology folder and compares:
* The peak flow 
* Time of peak flow

Inputs:
* Working model folder with pXX.hdf files and Hydrology folder with .dss files. 

Outputs: 
* Exports an HTML with interactive graphs of all ref line flow and junction flow of the model (should be easy to navigate between them both) X-axis time stamp, Y-axis Flow (CFS).

* Exports Excel table per HUC10 with the following columns:
    * Event Name from the .p text file
    * Plan number from the .hdf filename
    * Junction name (should match REF name as well) derived from the reference line first from the p##.hdf
    * Peak Flow (experienced ref line)
    * Time of Peak Flow (experienced ref line)
    * Peak Flow (expected junction flow)
    * Time of Peak Flow (expected junction flow)
    * Difference in Peak Flow (Reference - Junction)
    * Difference in Time of Peak Flow (Reference - Junction)

Folder Structure: 
* The output folder is created 2 folders above the input folder in "References" > [huc 10] > output files. It doesn't matter where the code folder lives, output folder is determined using input folder. 

Requirements: 
* regularhome.yml 