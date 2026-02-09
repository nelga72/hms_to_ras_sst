#This notebook contains all the user defined variables for the auto_bc_creation notebook. 
#It must be located a folder labelled src/inputs_files in the project directory structure.
import os
import pathlib as pl
import shutil

#set home directory - do not change
home = pl.Path(os.getcwd())


##############################################################
################  USER DEFINED VARIABLES  ####################
########### ONLY CHANGE THESE ACCORDING TO PROJECT   #########
##############################################################

model_prefix = 'wy_bh_' #This prefix is constant throughout each model in the project. Within the same project, the prefix should not be changing. 

#user to set variables for the project. tagged as parameter for papermill runs
project = 'wy_fy23'
huc = '1008001302' #this is the ds huc of the target huc10

schematic_input_type = '*.gdb'  #file type for the schematic input, either .gdb for geodatabases or .gpkg for geopackages.

temporary_overwrite = True #if True, you can limit the amount of target huc10s the code will be processed for. 
temporary_overwrite_us_huc = ['1008001301'] #list of the huc10s you would like to process. Only used if temporary_overwrite is True.

##############################################################
##############################################################


force = True #deletes outputs from any prior runs. Set to False to avoid this


#set up input and output folders
inputs = home/'inputs'
outputs_base = home/'outputs'
#create input and outputs folder if not already dirs
if not os.path.exists(inputs):
    os.makedirs(inputs)
if not os.path.exists(outputs_base):
    os.makedirs(outputs_base)

#hard variable for project
model_name = f'{model_prefix}{huc}' #outputs folder will be created with this name. This should be the same as the model name in the model file name.

#create output subfolder
if force:
    if os.path.exists(outputs_base/project/'eb_mod'/model_name):
        shutil.rmtree(outputs_base/project/'eb_mod'/model_name)

if not os.path.exists(outputs_base/project/'eb_mod'/model_name):
    os.makedirs(outputs_base/project/'eb_mod'/model_name)
    
outputs = outputs_base/project/'eb_mod'/model_name