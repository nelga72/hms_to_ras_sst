
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

#user to set variables for the project.
project = 'wy_fy23'  #project name, used to define input and output folders.
model_prefix = 'wy_bh_' #This prefix is constant throughout each model in the project. Within the same project, the prefix should not be changing. 
huc = '1008001303'  #HUC10 or model domain name for the model. This should be what is identified in the model file name. 
schematic_input_type = '*.gdb'  #file type for the schematic input, either .gdb for geodatabases or .gpkg for geopackages.

##############################################################
##############################################################

force = True #deletes outputs from any prior runs. Set to False to avoid this.

#define paths to input and output folders
inputs = home/'inputs'
outputs_base = home/'outputs'
#place all draft hydraulic models into the models_folder within inputs. Ensure that model files are located directly within the model's folder
models_folder = inputs/project/'hydraulic_models'

# Folder creation if not already present. Modifications not necessary but once this is imported, you can fill out folders with the necessary data.

#create required folders if not already present
if not os.path.exists(inputs):
    os.makedirs(inputs)
if not os.path.exists(models_folder):
    os.makedirs(models_folder)
#place all dictionaries in inputs/dictionaries folder, all files related to model events and the boundary condition go into the event_tie_ins folder
if not os.path.exists(inputs/project/'dictionaries'/'event_tie_ins'):
    os.makedirs(inputs/project/'dictionaries'/'event_tie_ins')
#place files provided from HMS related to hydrology into the inputs/project/hms_ras folder below. Ensure separated by model folder.
if not os.path.exists(inputs/project/'hms_ras'):   
    os.makedirs(inputs/project/'hms_ras')
#place files provided from HMS for schematics into the inputs/project/schematics folder below. Ensure the geodatabases are directly present
if not os.path.exists(inputs/project/'schematics'):
    os.makedirs(inputs/project/'schematics')

#outputs base folder
if not os.path.exists(outputs_base):
    os.makedirs(outputs_base)

#create output subfolder
model_name = f'{model_prefix}{huc}'
if force:
    if os.path.exists(outputs_base/project/model_name):
        shutil.rmtree(outputs_base/project/model_name)
if not os.path.exists(outputs_base/project/model_name):
    os.makedirs(outputs_base/project/model_name)
    
outputs = outputs_base/project/model_name

#create plot and notebook output location
#create output subfolder
if not os.path.exists(outputs_base/project/'plots'):
    os.makedirs(outputs_base/project/'plots')
plots = outputs_base/project/'plots'





