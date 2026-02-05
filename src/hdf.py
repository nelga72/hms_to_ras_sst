#imports
import zipfile
import h5py, os
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import pathlib as pl
import math
import pyproj
import shapely
from shapely import geometry
from shapely.geometry import Point, LineString
from shapely.ops import unary_union
from scipy.spatial import cKDTree
import re
import shutil
from pydsstools.heclib.dss import HecDss
from pydsstools.core import TimeSeriesContainer,UNDEFINED
from pydsstools.core import PairedDataContainer
from pydsstools.heclib.dss.HecDss import Open
from datetime import datetime
from collections import Counter
import json
import glob
import fiona

def get_current_ras_files(prj_file:pl.Path):
    """_summary_

    Args:
        prj_file (pl.Path): _description_

    Returns:
        _type_: _description_
    """
    #get current plan name and associated flow and geometry files.
    with open(prj_file,'r',encoding="ISO-8859-1") as prj:
        content = prj.read()
    current_plan = re.search('Current Plan=p\d*',content).group(0).split('=')[-1]
    assert current_plan is not None, "no plan listed as default in RAS"
    plan_path = prj_file.parent/f'{prj_file.stem}.{current_plan}'

    #get current plan, flow, geometry file and make copies in the outputs folder
    with open(plan_path) as plan_file:
        content = plan_file.read()
    in_geom = re.search('Geom File=g\d*',content).group(0).split('=')[-1]
    in_flow = re.search('Flow File=u\d*',content).group(0).split('=')[-1]
    assert in_geom is not None and in_flow is not None, "no geometry and flow listed as default in RAS"
    in_geom_path = prj_file.parent/f'{prj_file.stem}.{in_geom}'
    in_flow_path = prj_file.parent/f'{prj_file.stem}.{in_flow}'
    return plan_path, in_geom_path, in_flow_path

def get_ext_bc_lines(main_huc_perim:gpd.GeoDataFrame,geom_file:pl.Path,bc_line_offset:int):
    
    #geo hdf
    hf_ds = h5py.File(str(geom_file),'r')

    #get domain names
    domain_ds = str(list(hf_ds['Geometry']['2D Flow Areas']['Attributes'])[0][0]).strip("'b\'")

    #get global projection for all geometric operations. Assume the same between model plans
    crs_prj = str(hf_ds.attrs['Projection']).strip("'b\'")
    
    #get external boundary condition names for upstream models
    ext_bcs = [str(x[0]).strip("'b\'") for x in list(hf_ds['Geometry']['Boundary Condition Lines']['Attributes'])]

    #get geometry for external boundary conditions
    ext_faces = [list(y) for y in list(hf_ds['Geometry']['Boundary Condition Lines']['External Faces'])]

    #only get outflow bc face lines
    ds_bcs = {}
    for i, ext_bc in enumerate(ext_bcs):
        if ext_bc.lower().find('outflow')>=0:
            ds_bcs[ext_bc] = {'faces':[z[1] for z in ext_faces if z[0] == i],
                                  'fp_s':[z[2] for z in ext_faces if z[0] == i],
                                  'fp_e':[z[3] for z in ext_faces if z[0] == i],
                                 'start':[z[4] for z in ext_faces if z[0] == i],
                                  'stop':[z[5] for z in ext_faces if z[0] == i]}
    assert bool(ds_bcs), f"no external boundary conditions found in {geom_file}. Check geometry and that naming convention includes 'outflow_to_...'"
    
    #get index, geometry, and wsel of faces for entire mesh
    facepoint_coords = list(hf_ds['Geometry']['2D Flow Areas'][domain_ds]['FacePoints Coordinate'])
    ds_bc_df = pd.DataFrame({'name':pd.Series(dtype='str'),'face':pd.Series(dtype='int'),'fp_start':pd.Series(dtype='int'),'fp_end':pd.Series(dtype='int')})

    for ds_ext_bc_name, ds_data  in ds_bcs.items():
        #connect faces with face points
        ds_bc_df_tmp = pd.DataFrame({'name':[ds_ext_bc_name]*len(ds_data['faces']),'face':ds_data['faces'],'fp_start':ds_data['fp_s'],'fp_end':ds_data['fp_e']})
        ds_bc_df = pd.concat([ds_bc_df,ds_bc_df_tmp])

    #get coordinate information
    ds_bc_df['geometry'] = ds_bc_df.apply(lambda x: LineString([Point(facepoint_coords[x.fp_start]),Point(facepoint_coords[x.fp_end])]),axis=1)

    #trim start and end points
    for ds_ext_bc_name  in ds_bcs.keys():                 
        df = ds_bc_df.loc[ds_bc_df['name'] == ds_ext_bc_name]
        start_face, mod_start_array = [df.iloc[0].face, df.iloc[0].geometry.interpolate(-ds_bcs[ds_ext_bc_name]['start'][0])]
        ds_bc_df.loc[ds_bc_df['face'] == start_face, 'geometry'] = LineString([mod_start_array,Point(df.iloc[0].geometry.coords[-1])])
        stop_face, mod_stop_array = [df.iloc[-1].face, df.iloc[-1].geometry.interpolate(ds_bcs[ds_ext_bc_name]['stop'][-1]-ds_bcs[ds_ext_bc_name]['stop'][-2])]
        ds_bc_df.loc[ds_bc_df['face'] == stop_face, 'geometry'] = LineString([Point(df.iloc[-1].geometry.coords[0]),mod_stop_array])

    #convert to geodataframe
    ds_bc_gdf = gpd.GeoDataFrame(ds_bc_df,crs=crs_prj)

    #duplicate lines so that cell face end points can be recreated
    bc_wsels_pnts  = ds_bc_gdf.copy()
    bc_wsels_pnts['flag']= True
    ds_bc_gdf['flag']= False
    
    bc_wsels_pnts = pd.concat([bc_wsels_pnts,ds_bc_gdf.copy()],sort=True,axis=0)
    bc_wsels_pnts.geometry = bc_wsels_pnts.apply(lambda x: Point(x.geometry.coords[0]) if x.flag is True else Point(x.geometry.coords[1]),axis = 1)

    #prep setup
    bc_wsels_pnts['face_idx'] = bc_wsels_pnts.index
    bc_wsels_pnts.reset_index(inplace=True,drop=True)
    bc_wsels_pnts['point_idx'] = bc_wsels_pnts.index
    bc_wsels_pnts.reset_index(inplace=True,drop=True)

    #start
    df_line = pd.DataFrame(bc_wsels_pnts.iloc[0]).T
    running_dist = 0
    df_line['distance'] = 0
    used_points = [0]
    other_points = bc_wsels_pnts.loc[bc_wsels_pnts.point_idx.isin(used_points) == False]
    df_line['flag'] = df_line['flag'].astype(bool)
    while other_points.empty is False:
        id_pnt = used_points[-1]
        if other_points.loc[other_points.face_idx == bc_wsels_pnts.loc[id_pnt].face_idx].empty is False:
            connected_point = other_points.loc[other_points.face_idx == bc_wsels_pnts.loc[id_pnt].face_idx].copy()
            nA = np.array(list((bc_wsels_pnts.loc[id_pnt].geometry.x, bc_wsels_pnts.loc[id_pnt].geometry.y)))
            nB = np.array(list(connected_point.geometry.apply(lambda x: (x.x, x.y))))
            btree = cKDTree(nB)
            dist, idx = btree.query(nA, k=1)
            running_dist += dist
            connected_point['distance'] = running_dist
            df_line = pd.concat([df_line,connected_point],sort=False,axis=0)
            used_points.append(connected_point.point_idx.iloc[0])
        else:
            nA = np.array(list((bc_wsels_pnts.loc[id_pnt].geometry.x, bc_wsels_pnts.loc[id_pnt].geometry.y)))
            nB = np.array(list(other_points.geometry.apply(lambda x: (x.x, x.y))))
            btree = cKDTree(nB)
            dist, idx = btree.query(nA, k=1)
            running_dist += dist
            pnt_nearest = pd.DataFrame(other_points.iloc[idx]).T
            pnt_nearest['flag'] = pnt_nearest['flag'].astype(bool)
            pnt_nearest['distance'] = running_dist
            df_line = pd.concat([df_line,pnt_nearest],sort=False,axis=0)
            used_points.append(pnt_nearest.point_idx.iloc[0])
        other_points = bc_wsels_pnts.loc[bc_wsels_pnts.point_idx.isin(used_points) == False]

    #reformatting
    df_line.reset_index(inplace=True,drop=True)
    df_line['nominal_distance'] = df_line['distance']/df_line['distance'].max()

    #create line
    bc_line = LineString(np.array(list(df_line.geometry.apply(lambda x: (x.x, x.y)))))
    bc_line_df = pd.DataFrame()
    bc_line_df['name'] = [domain_ds,]
    bc_line_df['huc'] = re.search('\d+',pl.Path(geom_file).stem).group()
    bc_line_gdf = gpd.GeoDataFrame(bc_line_df, geometry=[bc_line,],crs=crs_prj)

    ### convert from faces to an offset line
    #offset lines and flip to confirm orientation
    ds_bc_gdf_left = bc_line_gdf.copy()
    ds_bc_gdf_right = bc_line_gdf.copy()
    ds_bc_gdf_left['geometry'] = ds_bc_gdf_left.apply(lambda x: x.geometry.parallel_offset(bc_line_offset, 'left'), axis=1)
    ds_bc_gdf_right['geometry'] = ds_bc_gdf_right.apply(lambda x: x.geometry.parallel_offset(bc_line_offset, 'right'), axis=1)

    #add function to flip geometry
    def reverse_line(geom):
        def _reverse(x, y):
            return x[::-1], y[::-1]

        return shapely.ops.transform(_reverse, geom)

    #flip lines that push to the interior when offset from the right
    inf_faces_int = ds_bc_gdf_right.loc[ds_bc_gdf_right.within(main_huc_perim)]

    inf_faces_flipped = bc_line_gdf.copy()
    if not inf_faces_int.empty:
        inf_faces_flipped.loc[list(inf_faces_int.index),'geometry'] = inf_faces_flipped.loc[list(inf_faces_int.index),'geometry'].apply(lambda x: reverse_line(x))

    # reoffset to the right (should all be exterior now
    ds_bc_gdf = inf_faces_flipped.copy()
    
    bc_line_gdf_offset = ds_bc_gdf.copy()
    bc_line_gdf_offset['perim_geometry'] = bc_line_gdf_offset['geometry']
    bc_line_gdf_offset['geometry'] = ds_bc_gdf.apply(lambda x: x.geometry.parallel_offset(bc_line_offset, 'right'), axis=1)
    return bc_line_gdf_offset

def get_int_bc_lines(geom_file:pl.Path,inflow_huc:str):
    
    #geo hdf
    hf_ds = h5py.File(str(geom_file),'r')

    #get domain names
    domain_ds = str(list(hf_ds['Geometry']['2D Flow Areas']['Attributes'])[0][0]).strip("'b\'")

    #get global projection for all geometric operations. Assume the same between model plans
    crs_prj = str(hf_ds.attrs['Projection']).strip("'b\'")
    
    #get external boundary condition names for current models
    ext_bcs = [str(x[0]).strip("'b\'") for x in list(hf_ds['Geometry']['Boundary Condition Lines']['Attributes'])]

    #get geometry for external boundary conditions
    ext_parts = [list(y) for y in list(hf_ds['Geometry']['Boundary Condition Lines']['Polyline Parts'])]
    
    #get geometry for external boundary conditions
    ext_points = [list(y) for y in list(hf_ds['Geometry']['Boundary Condition Lines']['Polyline Points'])]
    

    #only get outflow bc face lines
    ds_bcs = {}
    start = 0
    for i, ext_bc in enumerate(ext_bcs):
        if ext_bc.lower().find(f'inflow_from_{inflow_huc}')>=0:
            print(f'adding out of scope inflow from {inflow_huc}')
            ds_bcs[ext_bc] = ext_points[hf_ds['Geometry']['Boundary Condition Lines']['Polyline Parts'][i][0]+start:hf_ds['Geometry']['Boundary Condition Lines']['Polyline Parts'][i][1]+start]
        else:
            pass
            #print(f'WARNING: no upstream boundary condition found in {geom_file} for out of scope upstream huc {inflow_huc}. Check geometry and that naming convention includes "inflow_from_[huc#]"')
        start+= hf_ds['Geometry']['Boundary Condition Lines']['Polyline Parts'][i][1]

    assert bool(ds_bcs), f"no external boundary conditions found in {geom_file}. Check geometry and that naming convention includes 'inflow_from_...'"
    
    for inflow_bc, l_coords in ds_bcs.items():
        #create line
        bc_line = LineString(l_coords)
        bc_line_df = pd.DataFrame()
        bc_line_df['name'] = [inflow_bc,]
        inf_huc = inflow_bc.split('_')[-1]
        bc_line_df['huc'] = inf_huc
        bc_line_gdf = gpd.GeoDataFrame(bc_line_df, geometry=[bc_line,],crs=crs_prj)
        if inflow_bc == list(ds_bcs.keys())[0]:
            bc_line_gdfs = bc_line_gdf
        else:
            bc_line_gdfs = pd.concat([bc_line_gdfs,bc_line_gdf])

    return bc_line_gdf


##############hdf functions not used##########################
# #geo hdf
# hf_geo = h5py.File(str(geo_out_hdf),'a')

# #get domain name
# domain_geo = str(list(hf_geo['Geometry']['2D Flow Areas']['Attributes'])[0][0]).strip("'b\'")

# #get global projection for all geometric operations. Assume the same between model plans
# crs_prj = str(hf_geo.attrs['Projection']).strip("'b\'")

# inflow_dict = inflow_bcs.to_dict()

# exist_index = len(list(hf_geo['Geometry']['Boundary Condition Lines']['Polyline Info']))
# exist_len = len(list(hf_geo['Geometry']['Boundary Condition Lines']['Polyline Points']))
# polyline_points = list(hf_geo['Geometry']['Boundary Condition Lines']['Polyline Points'])
# exist_attributes_samp = hf_geo['Geometry']['Boundary Condition Lines']['Attributes'][0]

# att = []
# polyline_info = []
# polyline_parts = []

# for o in range(exist_index):
#     att.append(hf_geo['Geometry']['Boundary Condition Lines']['Attributes'][o])
#     polyline_info.append(hf_geo['Geometry']['Boundary Condition Lines']['Polyline Info'][o])
#     polyline_parts.append(hf_geo['Geometry']['Boundary Condition Lines']['Polyline Parts'][o])

# for i in list(inflow_dict['name'].keys()):
#     index = i+exist_index
#     shape = inflow_dict['geometry'][i]
#     array_len = len(np.array(inflow_dict['geometry'][i].coords))
#     attributes = exist_attributes_samp
#     attributes[0] = inflow_dict['name'][i]
#     attributes[1] = domain_geo.encode('ascii')
#     attributes[2] = 'External'
#     att.append(attributes)
#     #attributes[3] = shape.length
#     polyline_info.append(np.array([exist_len,array_len, index, 1]))
#     polyline_parts.append(np.array([0,array_len]))
#     for point in list(shape.coords):
#         polyline_points.append(np.array(point))
#     exist_len +=array_len
# del hf_geo['Geometry']['Boundary Condition Lines']['Attributes']
# del hf_geo['Geometry']['Boundary Condition Lines']['Polyline Info']
# del hf_geo['Geometry']['Boundary Condition Lines']['Polyline Parts']
# del hf_geo['Geometry']['Boundary Condition Lines']['Polyline Points']

# hf_geo.create_dataset('/Geometry/Boundary Condition Lines/Attributes', data=np.array(att))
# hf_geo.create_dataset('/Geometry/Boundary Condition Lines/Polyline Info', data=np.array(polyline_info))
# hf_geo.create_dataset('/Geometry/Boundary Condition Lines/Polyline Parts', data=np.array(polyline_parts))
# hf_geo.create_dataset('/Geometry/Boundary Condition Lines/Polyline Points', data=np.array(polyline_points))



def write_updated_ext_bc_files(domain_geo,dss_files,ext_index,ext_bc_dict,int_index,int_bc_dict,int_r_index,int_r_bc_dict, src_index,src_bc_dict, inputs,outputs,project,huc,start_id,src_huc_hms,model_name_prefix,append=False, force=False):
    #set variables
    bc_info = {}
    bc_info['in_plan_path']= ext_bc_dict['in_plan_path']
    bc_info['in_flow_path']= ext_bc_dict['in_flow_path']
    plans = []
    model_name = f'{model_name_prefix}{huc}' #folder within inputs folder    
    #get us subbasins for flow data
    with open(inputs/project/'dictionaries'/'Junction_Subbasins.json') as src:
        j_connect_sub = json.load(src)
    #write copies of the data 
    if not os.path.exists(outputs/f'{model_name}.prj') or force:
        prj_out_name = shutil.copy(inputs/project/'hydraulic_models'/huc/f'{model_name}.prj',outputs/f'{model_name}.prj')
    else:
        prj_out_name = outputs/f'{model_name}.prj'
    #geo_out_name = outputs/f'{model_name}.g{geo_id}'
    #geo_out_hdf = outputs/f'{model_name}.g{geo_id}.hdf'
    #ext_bc_dict['geo_id'] = geo_id 
    if append:
        #temp hardcode - need to rewrite this to read from the existing plan files.
        geo_id = 51
    else:
        geo_id = start_id
    bc_info['geo_id'] = geo_id
    #loop through dss files
    for dss_file in dss_files:
        #update id for each plan
        dss_name = pl.Path(dss_file).name
        dss_stem = pl.Path(dss_file).stem
        bc_info['dss_name'] = dss_name
        # set plan and flow names
        bc_info['flow_name'] = f'{dss_stem}_{start_id}'
        bc_info['plan_name'] = f'{dss_stem}_{start_id}'
        bc_info['out_num'] = start_id
        #copy into local model hydrology folder
        #if not os.path.exists(inputs/project/huc[4:8]/huc/'hydrology'/dss_name) or force:
        #    shutil.copy(dss_file,inputs/project/huc[4:8]/huc/'hydrology'/dss_name)

        fid = HecDss.Open(str(dss_file))
        pathname_pattern ="/*/*/*/*/*/*/"
        dss_list = fid.getPathnameList(pathname_pattern,sort=1)

        # print('currently at beginning of combining dss files',dss_list)

        #get relevant dss paths for external boundaries
        ext_bc_dict['dss_path'] = {i: None for i in ext_index}
        for i in ext_index:
            #get all relevant dss paths and condense using wildcard
            ext_junc = ext_bc_dict['junction'][i]
            ext_junc_finder = f'//{ext_junc}/FLOW-COMBINE/[\d\w\S]+'
            r = re.compile(ext_junc_finder)
            dss_matches = list(filter(r.match, dss_list))
            # print('Quick check for external boundary condition. remove print statement after testing:',dss_matches,src_huc_hms[ext_junc],huc)
            # ###################################################### MAJOR CHANGES REQUIRED. REVIEW.


            # ##added in function to copy upstream huc8 dss paths to the current dss file
            # if not dss_matches and src_huc_hms[ext_junc] != huc:
            #     print(f'us_huc8_triggered_external_junction {src_huc_hms[ext_junc]}_{ext_junc}')
            #     #dss_name_us = dss_name.replace('R','R-')
            #     dss_name_us = dss_name
            #     # assert os.path.exists(str(pl.Path(inputs/project/'us_dss'/f'HUC{huc[:8]}'/huc/f'us_dss_HUC{src_huc_hms[ext_junc][:8]}'/dss_name_us))), f'upstream dss file {dss_name_us} does not exist. Please request it from Dewberry for upstream huc8 {src_huc_hms[ext_junc][:8]}.'
            #     # fid_us_huc = HecDss.Open(str(pl.Path(inputs/project/'us_dss'/f'HUC{huc[:8]}'/huc/f'us_dss_HUC{src_huc_hms[ext_junc][:8]}'/dss_name_us)))
            #     assert os.path.exists(str(pl.Path(inputs/project/'transfer_dss'/f'{src_huc_hms[ext_junc][:8]}'/dss_name_us))), f'upstream dss file {dss_name_us} does not exist. Please request it from Dewberry for upstream huc8 {src_huc_hms[ext_junc][:8]}.'
            #     fid_us_huc = HecDss.Open(str(pl.Path(inputs/project/'transfer_dss'/f'{src_huc_hms[ext_junc][:8]}'/dss_name_us)))
            #     dss_list_us_huc = fid_us_huc.getPathnameList(pathname_pattern,sort=1)
            #     r = re.compile(ext_junc_finder)
            #     dss_matches = list(filter(r.match, dss_list_us_huc))
            #     #copy_path             
            #     for dss_match_us_huc in dss_matches:
            #         source_path = dss_match_us_huc
            #         target_path = dss_match_us_huc
            #         fid_us_huc.copyRecordsTo(fid,source_path,target_path)
            
            #print(ext_junc)
            assert len(dss_matches) > 0, f'cannot find a match for {ext_junc_finder} in {dss_name}'
            all_data = dss_matches[0]
            hg_all = all_data
            diff = [w for w in dss_matches[0].split('/') if w not in dss_matches[1].split('/')]
            for part in diff:
                hg_all = hg_all.replace(part,'*')
            ext_bc_dict['dss_path'][i] = hg_all
        #get relevant dss paths for internal boundaries
        int_bc_dict['dss_path'] = {i: None for i in int_index}
        # print('preparing the internal boundary conditions')
        # print('internal index',int_index)
        for i in int_index:
            # print(f'processing internal boundary condition {i}')
            int_junc = int_bc_dict['junction'][i]
            #get relevant subbasins to re-write flow values
            int_sbs = j_connect_sub[int_junc]
            # print('subbasin and junction',int_sbs,int_junc)
            #only add flow location information if there are contibuting subbasins
            # print('checking if internal junction has subbasins connected:',int_sbs)
            if int_sbs:
                int_junc_finder = f'//{int_junc}/FLOW-COMBINE/[\d\w\S]+'
                #get all relevant dss paths and condense using wildcard
                r = re.compile(int_junc_finder)
                dss_matches = list(filter(r.match, dss_list))
                # print('DSS matches found:', dss_matches)
                # ##NOT IMPLEMENTED function to copy upstream huc8 dss paths to the current dss file
                # if not dss_matches and src_huc_hms[int_junc] != huc:
                #     print('us_huc8_triggered_internal_junction')
                #     dss_name_us = dss_name.replace('R','R-')
                #     fid_us_huc = HecDss.Open(str(pl.Path(inputs/project/'us_dss'/f'HUC{huc[:8]}'/huc/f'us_dss_HUC{src_huc_hms[int_junc][:8]}'/dss_name_us)))
                #     dss_list_us_huc = fid_us_huc.getPathnameList(pathname_pattern,sort=1)
                #     r = re.compile(int_junc_finder)
                #     dss_matches = list(filter(r.match, dss_list_us_huc))
                #     #copy_path
                #     source_dss = fid_us_huc               
                #     target_dss = Open(str(inputs/project/huc[4:8]/huc/'hydrology'/dss_name))
                #     for dss_match_us_huc in dss_matches:
                #         source_path = dss_match_us_huc
                #         target_path = dss_match_us_huc
                #         source_dss.copyRecordsTo(target_dss,source_path,target_path)
                #     ##add in copy for dependent subbasins if needed
                #     ##
                #     ##
                
                
                assert len(dss_matches) > 0, f'cannot find a match for {int_junc_finder} in {dss_name}. Curtis needs to add in upstream HUC8 code. Ask him to do this.'
                all_data = dss_matches[0]
                hg_all = all_data
                diff = [w for w in dss_matches[0].split('/') if w not in dss_matches[1].split('/')]
                for part in diff:
                    hg_all = hg_all.replace(part,'*')
                #create new path name
                hg_name = hg_all[:hg_all.find('/FLOW-COMBINE/')]+'_SBS'+hg_all[hg_all.find('/FLOW-COMBINE/'):]
                int_bc_dict['dss_path'][i] = hg_name
                # print(hg_name)
                # print('create_junc_from_subs start')
                create_junc_from_subs(int_junc,int_sbs,fid, dss_list,hg_name)
                # print('create_junc_from_subs end')
        # print('int_index looped')
        
        #get relevant dss paths for reach internal boundaries
        int_r_bc_dict['dss_path'] = {i: None for i in int_r_index}
        # print('int_r_index',int_r_index)
        for i in int_r_index:
            # print(i)
            #get relevant junctions to re-write flow values
            reach = int_r_bc_dict['name'][i]
            us_junc = int_r_bc_dict['upstream'][i]
            subtraction = int(int_r_bc_dict['subtraction'][i])
            reduction = float(int_r_bc_dict['reduction'][i])

            # print(us_junc)

            #create a negative flow add flow location information if there are contibuting subbasins
            assert us_junc, f"no upstream junction connected to reach {int_r_bc_dict['name'][i]} with losses"
            r_junc_finder = f'//{us_junc}/FLOW-COMBINE/[\d\w\S]+'
            #get all relevant dss paths and condense using wildcard
            r = re.compile(r_junc_finder)
            dss_matches = list(filter(r.match, dss_list))
            all_data = dss_matches[0]
            hg_all = all_data
            # print(hg_all)
            diff = [w for w in dss_matches[0].split('/') if w not in dss_matches[1].split('/')]
            for part in diff:
                hg_all = hg_all.replace(part,'*')
            #create new path name
            hg_name = hg_all[:hg_all.find('/FLOW-COMBINE/')]+'_flow_reducer'+hg_all[hg_all.find('/FLOW-COMBINE/'):]
            int_r_bc_dict['dss_path'][i] = hg_name
            # print('made it this far', i)
            create_junc_from_reach(us_junc,dss_matches,subtraction,reduction,fid,dss_list,hg_name)
        # print('int_r_index looped')
        #get relevant dss paths for external boundaries
        src_bc_dict['dss_path'] = {i: None for i in src_index}
        for i in src_index:
            #get all relevant dss paths and condense using wildcard
            src_junc = src_bc_dict['junction'][i]
            src_junc_finder = f'//{src_junc}/FLOW/[\d\w\S]+'
            r = re.compile(src_junc_finder)
            dss_matches = list(filter(r.match, dss_list))
            all_data = dss_matches[0]
            hg_all = all_data
            diff = [w for w in dss_matches[0].split('/') if w not in dss_matches[1].split('/')]
            for part in diff:
                hg_all = hg_all.replace(part,'*')
            src_bc_dict['dss_path'][i] = hg_all
        # print('src_index looped')
        #get start and end times
        int_junc = int_bc_dict['junction'][0]
        int_junc_finder = f'//{int_junc}/FLOW-COMBINE/[\d\w\S]+'
        r = re.compile(int_junc_finder)
        dss_matches = list(filter(r.match, dss_list))
        start_date, end_date = get_ras_window_from_dss(fid,dss_matches)
        fid.close()
        
        #set plan file names
        plan_out_name = outputs/f'{model_name}.p{start_id}'
        # print(plan_out_name)
        flow_out_name = outputs/f'{model_name}.u{start_id}'
        #write updated flow file
        # print('writing updated flow file and plan file for plan id:',start_id)


        write_flow_file(domain_geo,bc_info,ext_bc_dict,ext_index,ext_bc_dict['in_flow_path'],flow_out_name,keep_exist_flow=False) #inflow path to pull first data from input folder and then the latter ones append to it.
        write_flow_file(domain_geo,bc_info,int_bc_dict,int_index,flow_out_name,flow_out_name)
        write_flow_file(domain_geo,bc_info,int_r_bc_dict,int_r_index,flow_out_name,flow_out_name)
        write_flow_file(domain_geo,bc_info,src_bc_dict,src_index,flow_out_name,flow_out_name)
        write_plan_file(domain_geo,bc_info,bc_info['in_plan_path'],plan_out_name, start_date,end_date) 
        plans.append(start_id)
        start_id +=1
    #write update project file
    with open(prj_out_name, "r+") as f:
        file_contents = f.read()
    if append:
        g_finder = 'Geom File=[\d\w]+\\n(?!Geom)'
        prj_flow_finder = 'Unsteady File=[\d\w]+\\n(?!Unsteady)'
        prj_plan_finder = 'Plan File=[\d\w]+\\n(?!Plan)'
    else:
        cp_finder = 'Current Plan=[\d\w]+\\n'
        cp_replace = re.search(cp_finder,file_contents)
        g_finder = 'Geom File=[\d\w]+\\n'
        prj_flow_finder = 'Unsteady File=[\d\w]+\\n' ##not needed but still included so that later code can run
        prj_plan_finder = 'Plan File=[\d\w]+\\n(?!Plan)'
        
    g_replace = re.search(g_finder,file_contents)
    prj_fl_replace = re.search(prj_flow_finder,file_contents)
    prj_pl_replace = re.search(prj_plan_finder,file_contents)
    # plan file updates
    if append:
        if g_replace.group().find(str(geo_id))>=0:
            new_prj = file_contents[:g_replace.span()[1]]+f'Geom File=g{geo_id}\n'+file_contents[g_replace.span()[1]:prj_fl_replace.span()[1]]
        else:
            new_prj = file_contents[:prj_fl_replace.span()[1]]
    else:
        new_prj = file_contents[:cp_replace.span()[0]]+f'Current Plan=p{plans[0]}\n'+file_contents[cp_replace.span()[1]:g_replace.span()[0]]+f'Geom File=g{geo_id}\n'
    for plan_id in plans:
        new_prj +=f'Unsteady File=u{plan_id}\n'
    if append:
        new_prj +=file_contents[prj_fl_replace.span()[1]:prj_pl_replace.span()[1]]
    for plan_id in plans:
        new_prj +=f'Plan File=p{plan_id}\n'
    new_prj += file_contents[prj_pl_replace.span()[1]:]
    with open(prj_out_name,'w') as p_o:
        p_o.write(new_prj)

    print('completion of entire function')


def write_plan_file(domain_geo,bc_info,in_plan_path,out_plan_path, start_date,end_date):
    plan_name = bc_info['plan_name'] 
    out_num = bc_info['out_num']
    geo_id = bc_info['geo_id']
    with open(in_plan_path, "r+") as f:
        file_contents = f.read()
    plan_flow_finder = 'Flow File=[\d\w]+'
    p_replace = re.search(plan_flow_finder,file_contents)
    # plan file updates
    sim_date = ','.join([start_date.split(' ')[0]]+[''.join(start_date.split(' ')[1].split(':')[:2])]+[end_date.split(' ')[0]]+[''.join(end_date.split(' ')[1].split(':')[:2])])
    new_plan = [f'Plan Title={plan_name}',
                f'Short Identifier={plan_name}',
                f'Simulation Date={sim_date}',
                f'Geom File=g{geo_id}',
                f'Flow File=u{out_num}']
    new_plan_txt = '\n'.join(new_plan)
    plan_out = new_plan_txt+file_contents[p_replace.span()[1]:]
    with open(out_plan_path,'w') as p_o:
        p_o.write(plan_out) 

def write_flow_file(domain_geo,bc_info,bc_dict,index,in_flow_path,out_flow_path,keep_exist_flow=True):
    dss_name = bc_info['dss_name'] 
    #read old text
    with open(in_flow_path, "r+") as f:
        file_contents = f.read()
    #write new text
    new_text = ''
    for i in index:
        junc = bc_dict['junction'][i]
        slope = bc_dict['slope'][i]
        path = bc_dict['dss_path'][i]
        ext_bc_template = [f'Boundary Location=                ,                ,        ,        ,                ,{domain_geo},                ,{junc},                                ',                                
            'Interval=1HOUR',
            'Flow Hydrograph= 0 ',
            'Stage Hydrograph TW Check=0',
            f'Flow Hydrograph Slope= {slope} ',
            f'DSS File=.\hydrology\{dss_name}',
            f'DSS Path={path}',
            'Use DSS=True',
            'Use Fixed Start Time=False',
            'Fixed Start Date/Time=,',
            'Is Critical Boundary=False',
            'Critical Boundary Flow=']
        new_text += '\n'.join(ext_bc_template)
        new_text += '\n'

    #find insertion points
    flow_title_finder= 'Flow Title=[\d\w]+'
    flow_start_finder =  'Friction Slope=[\d.,]+\\n(?![\d\D.,\n]+Friction Slope)'
    precip_start_finder = 'Met Point Raster Parameters=[\w\d.,]+\\n'
    precip_mode_finder = 'Precipitation Mode=\w+(?=\\n)'
    title_match = re.search(flow_title_finder,file_contents)
    flow_start_match = re.search(flow_start_finder,file_contents)
    precip_match = re.search(precip_start_finder,file_contents)
    precip_mode_match = re.search(precip_mode_finder,file_contents)
    flow_name = bc_info['flow_name']
    title = f'Flow Title={flow_name}'
    start = file_contents[title_match.span()[1]:flow_start_match.span()[1]:]
    exist_flow = file_contents[flow_start_match.span()[1]:precip_match.span()[0]:]
    precip_1 = file_contents[precip_match.span()[0]:precip_mode_match.span()[0]]
    precip_2 = 'Precipitation Mode=Disable'
    precip_3 = file_contents[precip_mode_match.span()[1]:]
    if keep_exist_flow:
        flow_out = title+start+exist_flow+new_text+precip_1+precip_2+precip_3
    else:
        flow_out = title+start+new_text+precip_1+precip_2+precip_3
    with open(out_flow_path, "w") as f_o:
        f_o.write(flow_out)

def update_flow_file(dss_name,dss_path,in_flow_path,outflow_huc,domain_name):
    with open(in_flow_path, "r") as f:
        file_contents = f.read()
    #edit if flow file 
    ext_bc_template_update = [f'Interval=1HOUR',
        'Rating Curve= 0 ',
        f'Flow Hydrograph Slope= 0 ',
        f'DSS File=.\Hydrology\{dss_name}',
        f'DSS Path={dss_path}',
        'Use DSS=True',
        'Use Fixed Start Time=False',
        'Fixed Start Date/Time=,',
        'Is Critical Boundary=False',
        'Critical Boundary Flow=']

    #find insertion points
    flow_update_start =  domain_name+'[,\s]*[Oo]utflow[ _]to[ _]'+outflow_huc+'[\d\s\w\S]*?\n'
    flow_update_end = domain_name+'[,\s]*[Oo]utflow[ _]to[ _]'+outflow_huc+'[\d\s\w\S]*?\n(?=[Boundary Location])'
    flow_start_match = re.search(flow_update_start,file_contents)
    flow_end_match = re.search(flow_update_end,file_contents)
    assert flow_update_start, 'check naming convention of outflow bc'
    flow_update_out = file_contents[:flow_start_match.span()[1]]+'\n'.join(ext_bc_template_update)+'\n'+file_contents[flow_end_match.span()[1]:]
    with open(in_flow_path, "w") as f_o:
        f_o.write(flow_update_out)

def update_flow_file_stage(dss_name,dss_path,in_flow_path,outflow_huc,domain_name):
    with open(in_flow_path, "r") as f:
        file_contents = f.read()
    #edit if flow file 
    ext_bc_template_update = [f'Interval=1HOUR',
        'Stage Hydrograph= 0 ',
        'Stage Hydrograph Use Initial Stage=-1',
        'Stage Hydrograph TW Check=0',
        f'Flow Hydrograph Slope= 0 ',
        f'DSS File=.\Hydrology\{dss_name}',
        f'DSS Path={dss_path}',
        'Use DSS=True',
        'Use Fixed Start Time=False',
        'Fixed Start Date/Time=,',
        'Is Critical Boundary=False',
        'Critical Boundary Flow=']

    #find insertion points
    flow_update_start =  domain_name+'[,\s]*[Oo]utflow[ _]to[ _]'+outflow_huc+'[\d\s\w\S]*?\n'
    flow_update_end = domain_name+'[,\s]*[Oo]utflow[ _]to[ _]'+outflow_huc+'[\d\s\w\S]*?\n(?=[Boundary Location])'
    flow_start_match = re.search(flow_update_start,file_contents)
    flow_end_match = re.search(flow_update_end,file_contents)
    assert flow_update_start, 'check naming convention of outflow bc'
    flow_update_out = file_contents[:flow_start_match.span()[1]]+'\n'.join(ext_bc_template_update)+'\n'+file_contents[flow_end_match.span()[1]:]
    with open(in_flow_path, "w") as f_o:
        f_o.write(flow_update_out)

def create_bc_from_junc(bc_connections,fid,fid_us,stage_data_ds,flow_data_ds,dss_matches_us,rating_or_stage,event_us_flow_expected):
    """Writes a rating curve boundary condition to a dss file after reading information from a stage/flow hydrograph from a separate dss file.
    Parameters
    ----------
    bc_connections :
        a dictionary providing the huc, junctions, and dss path related to a specific recurrence interval
    fid :
        open pydsstool pointing to the downstream dss output file that will be read
    fid_us :
        open pydsstool pointing to the upstream dss output file that will be written to
    stage_data_ds :
        a list of the stage data dss paths that will be read from for the downstream dss file
    flow_data_ds :
        a list of the flow data dss paths that will be read from for the downstream dss file
    dss_matches_us :
        A list of dss matches from the upstream event dss (fid_us)
    rating_or_stage :
        A dictionary with the model names and their respective downstream boundary condition
    event_us_flow_expected :
        The expected flow for the selected event for the current specific recurrence interval
    Returns
    -------
    String
        A path name to the dss file path where the new rating curve information was written to.
    """
    hg_all = stage_data_ds[0]
    diff = [w for w in stage_data_ds[0].split('/') if w not in stage_data_ds[1].split('/')]
    for part in diff:
        hg_all = hg_all.replace(part,'*')
    #read data from ds huc 
    diffs={}
    #get date differences  
    for dss_match in stage_data_ds:
        for compare in stage_data_ds:
            for part in dss_match.split('/'):
                if part not in compare.split('/'):
                    i = dss_match.split('/').index(part)
                    if i not in diffs.keys():
                        diffs[i] = [part]
                    else:
                        if part not in diffs[i]:
                            diffs[i] = diffs[i]+[part]
   
    #get stage time series information
    stage_series = []
    for time_w in diffs[4]:    #hardcoded in part 4 of diffs to deal with time window. Need to assess whether there is a better way to do this
        bc_finder = f'[\d\w\S\s]+{time_w}[\d\w\S]+'
        r = re.compile(bc_finder)
        dss_match = list(filter(r.match, stage_data_ds))
        assert len(dss_match) != 0, f"there is no dss path matches for {r}{stage_data_ds}"
        assert len(dss_match) == 1, f"there is more than one dss path with the same bc junction and time window. debug {dss_match}"
        tss = fid.read_ts(dss_match[0])
        f = np.array(tss.values, dtype=float)
        assert None not in f, f'there are Nones in {dss_match[0]}'
        f[f<0] = 0
        if f.size == 1:
            if not np.isnan(f):
                stage_series.append(f)
            else:
                print(f'{dss_match} does not have data')
        else:
            stage_series.append(f)
    assert stage_series, f'data is not available for {stage_data_ds}'
    #find non null range of stage values and save to trim flow values as well
    stages_list = np.concatenate(stage_series).ravel()
    stages_non_null = np.where(stages_list != 0)[0]
    start_non_null = stages_non_null[0]
    stop_non_null = stages_non_null[-1]
    stages_null = np.where(stages_list == 0)[0]
    stages = stages_list[start_non_null:stop_non_null]
    
    #get flow time series
    flow_series = []
    for time_w in diffs[4]:   #hardcoded in part 4 of diffs to deal with time window. Need to assess whether there is a better way to do this
        bc_finder = f'[\d\w\S\s]+{time_w}[\d\w\S]+'
        r = re.compile(bc_finder)
        dss_match_f = list(filter(r.match, flow_data_ds))
        assert len(dss_match_f) != 0, f"there is no dss path matches for {r}{flow_data_ds}"
        assert len(dss_match_f) == 1, "there is more than one dss path with the same bc junction and time window. debug"
        tsf = fid.read_ts(dss_match_f[0])
        f = np.array(tsf.values, dtype=float)
        #f = np.copy(tsf.values)
        f[f<0] = 0
        if f.size == 1:
            if not np.isnan(f):
                flow_series.append(f)
            else:
                pass
        else:
            flow_series.append(f)
    assert flow_series, f'data is not available for {flow_data_ds}'
    flows_list = np.concatenate(flow_series).ravel()
    flows =  flows_list[start_non_null:stop_non_null]
    #confirm that clipping of np array didn't cause any inconsistencies
    assert len(stages) == len(flows), 'trimming of dss data resulted in mismatch stage / flow indices. Talk to curtis to update the code'    
    #analysis to create rating curve based on trial and error .
    #return stages, flows
    ##get rising and falling limbs of hydrograph without any looping
    stages_incline = []
    flows_incline = []
    stages_decline = []
    flows_decline = []
    min_flow = flows.min()
    max_stage = stages.max() #get max elevation
    max_stage_index = np.where(stages == stages.max())[0][0]
    max_flow_c = flows[max_stage_index] #get associated flow value
    # #filter hydrograph to only include flow and stage changes
    # change = np.intersect1d(np.where(np.gradient(stages) != 0),np.where(np.gradient(flows) != 0))
    # sub_flows = flows[change]
    # sub_stages = stages[change]

    #no_change = len(np.intersect1d(np.where(np.gradient(stages) == 0),np.where(np.gradient(flows) == 0))) 
    # #remove counts less than 5 cfs but falling stage hydrograph
    # zero_flow = np.where(flows<5)
    # stage_grad = np.gradient(stages) != 0
    # zero_flow_removals = len(np.where(stage_grad[zero_flow] == True)[0])

    #pull data that would be used for rating curve
    last_stage = stages[0]
    last_flow = flows[0]
    for stage, flow in zip(stages[:max_stage_index+1],flows[:max_stage_index+1]):
        if stage > last_stage and flow > last_flow:
            stages_incline.append(stage)
            flows_incline.append(flow)
            last_stage, last_flow = stage, flow
    last_stage = max_stage
    last_flow = max_flow_c
    for stage, flow in zip(stages[max_stage_index:],flows[max_stage_index:]):
        if stage < last_stage and flow < last_flow:
            stages_decline.append(stage)
            flows_decline.append(flow)
            last_stage, last_flow = stage, flow

    #ideally we want the falling limb unless the plan doesn't capture it well
    #rule 1 if one of the limbs doesn't exist, use the other one. assert that at least one exists

    if len(flows_decline) + len(flows_incline) == 0:
        print('no rising or falling data available, ponding condition will be triggered')
        side = 'none'
        stages_falling = np.array([])
        flows_falling = np.array([])
    else:
        if len(flows_incline) == 0:
            stages_falling = np.flip(np.array(stages_decline))
            flows_falling = np.flip(np.array(flows_decline))
            side='falling'
        elif len(flows_decline) == 0:
            stages_falling = np.array(stages_incline)
            flows_falling = np.array(flows_incline)
            side='rising'
        else:
            #rule 2 if min flow on falling limb is within  100 cfs of the min flow , then use it no matter what
            if np.array(flows_decline).min() <= (100+min_flow):
                stages_falling = np.flip(np.array(stages_decline))
                flows_falling = np.flip(np.array(flows_decline))
                side='falling'
            else:
                #overrule the rising limb if the number of data points is much greater on the falling limb.
                #this signals that while the rising limb has a lower start flow, the falling limb still has enough data to use for a rating curve.
                if len(stages_decline)*0.05 > len(stages_incline):
                        stages_falling = np.flip(np.array(stages_decline))
                        flows_falling = np.flip(np.array(flows_decline))
                        side='falling'
                else:
                    #rising limb has min elevation and enough points to be used for a rating curve
                    stages_falling = np.array(stages_incline)
                    flows_falling = np.array(flows_incline)
                    side='rising'
    
    if side == 'rising':
        stages_sub = stages[:max_stage_index+1]
        flows_sub = flows[:max_stage_index+1]
        flows_sub = flows_sub[np.where(flows_sub >= 1)]
        stages_sub = stages_sub[np.where(flows_sub >= 1)]
    elif side == 'falling':
        stages_sub = stages[max_stage_index:]
        flows_sub = flows[max_stage_index:]
        flows_sub = flows_sub[np.where(flows_sub >= 1)]
        stages_sub = stages_sub[np.where(flows_sub >= 1)]
    else:
        stages_sub = np.array([])
        flows_sub = np.array([])
    
    ##removed in place of pre-determined approach from dictionary but could be refined and used for future projects
    # #check for backwater ponding from downstream lakes. This check flags for backwater control if more than 5% of the flow increases / decreases do not result in corresponding stage changes
    # if len(stages_falling) < 3 or len(stages_sub) < 3:
    #     gradient_thresh = 0
    #     print('WARNING: ponding condition triggered for', dss_matches_us[0], 'because of limited rating curve data')
    # else:
    #     #create weighted gradient threshold parameter to identify where flow and stage gradients don't align, with higher priority given to higher stages.
    #     depth_stages_sub = np.subtract(stages_sub,stages_sub.min())
    #     norm_stages_sub = np.divide(depth_stages_sub,depth_stages_sub.max())
    #     gradient_thresh = (np.multiply(np.equal(np.gradient(stages_sub) > 0, np.gradient(flows_sub) > 0),norm_stages_sub).sum())/norm_stages_sub.sum()
    #     if gradient_thresh < 0.75:
    #         #debug
    #         out_path = pl.Path(os.getcwd())/'outputs'/'plots'
    #         if not os.path.exists(out_path):
    #             os.makedirs(out_path)
    #         fig, axs = plt.subplots(3, 1)
    #         axs[0].plot(stages_sub)
    #         axs[0].set_title('stages')
    #         axs[1].plot(flows_sub)
    #         axs[1].set_title('flows')
    #         axs[2].plot(flows_falling,stages_falling)
    #         axs[2].set_title('Rating Curve')
    #         name = dss_matches_us[0].split('/')[2]+'_'+dss_matches_us[0].split('/')[-2].replace(':','-')
    #         fig.savefig(out_path/f'{name}.png')
    #         plt.close()
    #         #np.save(f'C:\\temp\\{name}_flow.npy',flows)
    #         #np.save(f'C:\\temp\\{name}_stages.npy',stages)
    #         print('WARNING: ponding condition triggered for', dss_matches_us[0], 'because of lack of stage discharge relationship')
    #         print(gradient_thresh)
    # #print(dss_matches_us[0])
    # if gradient_thresh < 0.75:

    if rating_or_stage == "stage_hydrograph":
        #assume ponding control and create stage hydrograph boundary condition instead
        #print warning if  the peak flow is close to the peak stage not controlling the peak stage. Assume 4 hours timestep difference at least.
        #if abs(max_stage_index - np.where(flows == flows.max())[0][0]) > 1:
            #print('warning, double check that ponding is correctly flagged for this event') 

        #make name to indicate stage hydrograph instead of rating curve
        tsc = TimeSeriesContainer()
        tsc.interval = tss.interval
        tsc.units = tss.units
        tsc.type = tss.type
        
        #read temp data from us huc 
        diffs={}
        #get date differences  
        for dss_match in dss_matches_us:
            for compare in dss_matches_us:
                for part in dss_match.split('/'):
                    if part not in compare.split('/'):
                        i = dss_match.split('/').index(part)
                        if i not in diffs.keys():
                            diffs[i] = [part]
                        else:
                            if part not in diffs[i]:
                                diffs[i] = diffs[i]+[part]

        #create new path name for ts adjusted data
        hg_name = hg_all[:hg_all.find('/STAGE/')]+'_stage_hyd'+hg_all[hg_all.find('/STAGE/'):]
        
        hg_all_us = dss_matches_us[0]
        diff_us = [w for w in dss_matches_us[0].split('/') if w not in dss_matches_us[1].split('/')]
        for part in diff_us:
            hg_name_us = hg_all_us.replace(part,'*')
        
        tsfu = fid_us.read_ts(hg_name_us)
        #print(tsfu.startDateTime)
        f = np.array(tsfu.values, dtype=float)
        assert None not in f, f'there are Nones in {dss_matches_us[0]}'
        f[f<0] = 0
        if f.size == 1:
            if np.isnan(f):
                print(f'{dss_matches_us} does not have data')
                
        #break into month chunks
        f_data = []
        mon_len = [0]
        st_month = tsfu.pytimes[0].month
        mon = [st_month]
        count = 0
        for dt in tsfu.pytimes:
            if dt.month == st_month:
                count+=1
            else:
                mon_len.append(count)
                st_month = dt.month
                mon.append(st_month)
                f_data.append(f[mon_len[-2]:mon_len[-1]])
        f_data.append(f[mon_len[-1]:])
        offset = tsfu.pytimes[0].strftime("%d")
        start = 0
        end = 0
        if len(diffs[4]) != len(f_data):
            print(f"warning. a timestep is encroaching on a new month {diffs[4]} vs {tsfu.pytimes[-1]}. Automatically chopping off this timeslice ")
        for time_w, f in zip(diffs[4],f_data[:len(diffs[4])]):  
            #print(time_w)

            #assert f, f'data is not available for {stage_data_ds}'
            #find non null range of flow values and save to trim flow values as well
            flow_list = f.ravel()
            flow_non_null = np.where(flow_list != 0)[0]
            start_non_null = flow_non_null[0]
            stop_non_null = flow_non_null[-1]
            
            tsc.pathname = hg_name.replace('*',time_w)

            #pad start time forward based on flow data start. Only applies to first time window
            if start == 0:
                start_pad = start_non_null//(tsc.interval//tsfu.interval)
                #print(f'padding start of time series ({start_pad//24} days and {start_pad%24} hours)')
                stages = np.pad(stages,(int(start_pad),0))
                end+=start_pad
                #align to top of hour start day
                tsc.startDateTime = offset+time_w[2:] +' 01:00:00'
            else:
                #align to top of hour 1st day
                tsc.startDateTime = time_w +' 01:00:00'
            
            assert tsc.interval//tsfu.interval == 4, 'code assumes output mapping interval of 1 hour and flow input of 15 min'
            #end of slice based on end of flow time window
            end += int(len(flow_non_null)//(tsc.interval/tsfu.interval))+1
            #handle if slice isn't shorter than flow time window
            if len(stages[start:end]) < end-start:
                #print('padding end of ts')
                if stages[start:].any():
                    a = np.pad(stages[start:],(0,end-start-len(stages[start:end])+1),mode='edge')
                else:
                    a = np.pad(np.array([stages[-1]]),(0,end-start-len(stages[start:end])),mode='edge')
                end+=1
                a[a==0] = np.nan
                tsc.values = a
            else:
                a = stages[start:end]
                a[a==0] = np.nan
                tsc.values = a
            tsc.numberValues = end-start
            #return tsc
            #print(tsc.startDateTime)
            assert len(tsc.values) == tsc.numberValues, f"{end-start} and {len(tsc.values)} do not equal {tsc.numberValues} for {ri} {event}"
            fid_us.deletePathname(tsc.pathname)
            fid_us.put_ts(tsc)
            start=end
        return hg_name
    else:
        assert rating_or_stage == "rating_curve", f"{rating_or_stage} was not expected. Should be rating_curve or stage_hydrograph"
        #create new path name rating curve
        hg_parts = hg_all.split('/')
        hg_name = '/'.join([hg_parts[0],hg_parts[1],f'{hg_parts[2]}_rc','-','','',hg_parts[-2],hg_parts[-1]])
        
        
        #reset datapoints
        max_stage = stages_falling.max() #get max elevation
        max_stage_index = np.where(stages_falling == stages_falling.max())[0][0]
        max_flow_c = flows_falling[max_stage_index] #get associated flow value
        min_flow = flows_falling.min()
        rc_flow = [0]
        #use 30th percentile flow instead of using lowest stage to avoid bad low flow velocities
        temp_flow = np.copy(flows_falling)
        temp_flow[temp_flow == 0] = np.nan
        perc_flow = np.nanpercentile(temp_flow, 30,method='closest_observation')
        perc_flow_index = np.where(flows_falling == perc_flow)[0][0]
        #print(perc_flow,stages_falling[perc_flow_index],stages_falling.min())
        rc_stage_tmp = [stages_falling[perc_flow_index]]
        max_flow_offset_l = max_flow_c *0.2
        max_flow_offset_u = max_flow_c *1.8

        ## add datapoints for low flow
        inflection_points = [1,10,30,60,100,200, 400, 800, 2000]
        inflection_points_capped = [min_flow + num for num in inflection_points if min_flow + num < max_flow_offset_l]
        for value in inflection_points_capped:
            stage_interp = np.interp(value,flows_falling,stages_falling)
            if stage_interp > rc_stage_tmp[-1]:
                rc_flow.append(value)
                rc_stage_tmp.append(stage_interp) 

        #add max points + steeper last point for stability
        #stage
        rc_stage_tmp+=[max_stage-.5,max_stage,max_stage+1]
        rc_stage = np.array([rc_stage_tmp],dtype=np.float32) #np.array([[min_stage,max_stage-0.5,max_stage,max_stage+0.5]],dtype=np.float32)

        #flow
        rc_flow+=[max_flow_offset_l,max_flow_c,max_flow_offset_u]
        rc_flow_ratio = (event_us_flow_expected/max_flow_c)
        rc_flow = [item * rc_flow_ratio for item in rc_flow]
        #save to dss
        pdc = PairedDataContainer()
        pdc.pathname = hg_name
        pdc.curve_no = len(rc_stage)
        pdc.independent_axis = rc_flow
        pdc.data_no = len(rc_stage[0])
        pdc.curves = rc_stage #np.array([[7000,7001]],dtype=np.float32)
        pdc.labels_list = ['event']
        pdc.independent_units = '(cfs)'
        pdc.dependent_units = '(ft)'

        fid_us.put_pd(pdc)
        fid_us.close()
        return hg_name
   
        
        
def add_inflow_bcs(in_geom_path,index,bc_dict,domain_geo, geo_out,add_breaklines = False):
    #read old text
    with open(in_geom_path, "r+") as f:
        file_contents = f.read()
    new_text = ''
    breaklines = ''
    for i in index:
        #add bc lines (which are offset)
        geom = bc_dict['geometry'][i]
        first_coord = ' , '.join(str(coord) for coord in geom.coords[0])
        mid_coord = ' , '.join(str(coord) for coord in geom.line_interpolate_point(geom.length/2).coords[0])
        last_coord = ' , '.join(str(coord) for coord in geom.coords[-1])
        coord_len = len(geom.coords)
        junc = bc_dict['junction'][i]
        ext_bc_template = [f'BC Line Name={junc}',
                f'BC Line Storage Area= {domain_geo}',
                f'BC Line Start Position={first_coord}', 
                f'BC Line Middle Position={mid_coord}', 
                f'BC Line End Position={last_coord} ',
                f'BC Line Arc= {coord_len}']
        for i, coord in enumerate(geom.coords):
            if i % 2 == 0:
                coord_str = f"{coord[0]:16f}"+f"{coord[1]:16f}"
            else:
                coord_str += f"{coord[0]:16f}"+f"{coord[1]:16f}"
                ext_bc_template.append(coord_str)
                coord_str = ''
        if coord_str != '':
            ext_bc_template.append(coord_str)
        ext_bc_template.append('BC Line Text Position= 1.79769313486232E+308 , 1.79769313486232E+308')
        new_text += '\n'.join(ext_bc_template)
        new_text += '\n'
        #add breaklines
        if add_breaklines:
            bl_geom = bc_dict['perim_geometry'][i]
            bl_len = len(bl_geom.coords)
            bl_template = [f'BreakLine Name=bl_{junc}',
                    'BreakLine CellSize Min=',
                    'BreakLine CellSize Max=',
                    'BreakLine Near Repeats=0',
                    'BreakLine Protection Radius=0',
                    f'BreakLine Polyline= {bl_len}']
            for i, coord in enumerate(bl_geom.coords):
                if i % 2 == 0:
                    bl_coord_str = f"{coord[0]:16f}"+f"{coord[1]:16f}"
                else:
                    bl_coord_str += f"{coord[0]:16f}"+f"{coord[1]:16f}"
                    bl_template.append(bl_coord_str)
                    bl_coord_str = ''
            if bl_coord_str != '':
                bl_template.append(bl_coord_str)
            breaklines += '\n'.join(bl_template)
            breaklines += '\n'

    #find insertion points
    geo_bc_finder= 'LCMann Time'
    geo_bl_finder =  'BC Line Name'
    bc_insertion = re.search(geo_bc_finder,file_contents)
    bl_insertion = re.search(geo_bl_finder,file_contents)
    if not bl_insertion:
        bl_insertion = bc_insertion
    exist_bc = file_contents[bl_insertion.start():bc_insertion.start()]
    start = file_contents[:bl_insertion.start()]
    if add_breaklines:
        middle = breaklines+exist_bc+new_text
    else:
        middle = exist_bc+new_text
    end = file_contents[bc_insertion.start():]
    geo_out_str = start+middle+end

    with open(geo_out, "w") as g_o:
        g_o.write(geo_out_str)

def create_junc_from_subs(int_junc,int_sbs,fid,dss_list,hg_name):
    # print('starting variable input',int_sbs)
    diffs = {}
    #assume subbasins have the same time windows
    #get time windows
    int_sb_base = int_sbs[0]
    int_sb_finder = f'//{int_sb_base}/FLOW/[\d\w\S]+'
    # print(int_sb_finder)
    r = re.compile(int_sb_finder)
    dss_matches = list(filter(r.match, dss_list))
    # print('current dss matches',dss_matches)
    #get date differences  
    for dss_match in dss_matches:
        for compare in dss_matches:
            for part in dss_match.split('/'):
                if part not in compare.split('/'):
                    i = dss_match.split('/').index(part)
                    if i not in diffs.keys():
                        diffs[i] = [part]
                    else:
                        if part not in diffs[i]:
                            diffs[i] = diffs[i]+[part]
    #hardcoded in part 4 dealing with time window. Need to assess whether there is a better way to do this
    # print(diffs)
    for time_w in diffs[4]:
        flow_values = []
        for int_sb in int_sbs:
            # print(time_w) 
            int_sb_finder = f'//{int_sb}/FLOW/{time_w}/[\d\w\S]+'
            r = re.compile(int_sb_finder)
            dss_matches = list(filter(r.match, dss_list))
            assert len(dss_matches) == 1, "there is more than one dss path with the same subbasin and time window. debug"
            ts = fid.read_ts(dss_matches[0])
            f = np.copy(ts.values)
            f[f<0] = 0
            flow_values.append(f)
        temp = fid.read_ts(dss_matches[0])
        tsc = TimeSeriesContainer()
        tsc.pathname = hg_name.replace('*',time_w)
        tsc.startDateTime = temp.startDateTime
        tsc.numberValues = temp.numberValues
        tsc.units = temp.units
        tsc.type = temp.type
        tsc.interval = temp.interval
        if len(int_sbs) > 1:
            tsc.values = sum(flow_values)
            #tsc.nodata = ~np.array(sum(flow_values)).astype(bool)
        else:
            tsc.values = flow_values[0]
        fid.deletePathname(tsc.pathname)
        fid.put_ts(tsc)

def create_junc_from_reach(us_junc,dss_matches,subtraction,reduction,fid,dss_list,hg_name):
    #reach,us_junc,fid, dss_list,hg_name
    diffs = {}
    #get date differences   
    for dss_match in dss_matches:
        for compare in dss_matches:
            for part in dss_match.split('/'):
                if part not in compare.split('/'):
                    i = dss_match.split('/').index(part)
                    if i not in diffs.keys():
                        diffs[i] = [part]
                    else:
                        if part not in diffs[i]:
                            diffs[i] = diffs[i]+[part]
    #hardcoded in part 4 dealing with time window. Need to assess whether there is a better way to do this
    for time_w in diffs[4]:
        flow_values = []
        us_j_finder = f'//{us_junc}/FLOW-COMBINE/{time_w}/[\d\w\S]+'
        r = re.compile(us_j_finder)
        dss_matches = list(filter(r.match, dss_list))
        assert len(dss_matches) == 1, "there is more than one dss path with the same subbasin and time window. debug"
        ts = fid.read_ts(dss_matches[0])
        f = np.copy(ts.values)
        f[f<0] = 0
        flow_values.append(f)
        temp = fid.read_ts(dss_matches[0])
        tsc = TimeSeriesContainer()
        tsc.pathname = hg_name.replace('*',time_w)
        tsc.startDateTime = temp.startDateTime
        tsc.numberValues = temp.numberValues
        tsc.units = temp.units
        tsc.type = temp.type
        tsc.interval = temp.interval
        #math for reduction
        reduction_2 = np.subtract(f,np.full(f.shape, subtraction))
        reduction_2[reduction_2<0] = 0
        reduction_3 = np.multiply(reduction_2,np.full(f.shape, -(reduction)))
        tsc.values = reduction_3
        fid.deletePathname(tsc.pathname)
        fid.put_ts(tsc)
        
def get_ras_window_from_dss(fid,dss_paths):
    start_date = datetime.strptime('01Jan9999', '%d%b%Y') #fake start date way in the future
    end_date = datetime.strptime('01Jan0001', '%d%b%Y') #fake end date way in the past
    for path_data in dss_paths:
        ts = fid.read_ts(path_data,trim_missing=True)
        if datetime.strptime(ts.startDateTime[:9], '%d%b%Y') < start_date:
            start_date = datetime.strptime(ts.startDateTime[:9], '%d%b%Y')
            start_date_ras = ts.startDateTime
        if datetime.strptime(ts.endDateTime[:9], '%d%b%Y') > end_date:
            end_date = datetime.strptime(ts.endDateTime[:9], '%d%b%Y')
            end_date_ras = ts.endDateTime
    assert start_date != datetime.strptime('01Jan9999', '%d%b%Y') and  end_date != datetime.strptime('01Jan0001', '%d%b%Y'), "start and end dates are bad"
    return start_date_ras, end_date_ras

def get_sst_storms_by_recurrence(huc,hms_shps,flag=True):
    events_dict = {huc:{}}
    unique_events = []
    short_ids_dict = {}
    for ri_shp in hms_shps:
        f = re.search('_\d+.\d+[pm]*_',ri_shp)
        ri = f.group()[1:-1]
        events_dict[huc][ri] = []
        gdf = gpd.read_file(ri_shp)
        cols = list(gdf.columns.to_list())

        #################################################################
        #may require edits if the geojson being read changes, or if the naming convention for storms changes        
        r = re.compile('P\d+_R-Y\d+-E\d+')
        #################################################################
        col_matches = list(filter(r.match, cols))
        if col_matches:
            events_dict[huc][ri] = col_matches
            for col in col_matches:
                unique_events.append(col)
                r_short = re.search('Y\d+-E\d+',col)
                short_ids_dict[col] = r_short.group()
        else:
            events_dict[huc][ri] = []
    counts = Counter(unique_events)
    if counts.most_common()[0][1] <= 1:
        print(f"Note that the same storm event {counts.most_common()[0][0]} is used for more than recurrence interval")
    short_counts = Counter(list(short_ids_dict.values()))
    # if flag: 
    #     assert short_counts.most_common()[0][1] <= 1, f"the same storm event number and year number {short_counts.most_common()[0][0]} is used for more than one realization. Discuss with Matt D about adding the RX to the dss file names"
    return events_dict

def add_supplemental_dss():
    return None

# def get_sst_storms_by_recurrence_ds_huc(events_dict,huc,htmls,junction,j_to_j):
#     original_j = junction
#     for html in htmls:
#         html_file = open(html, "r") 
#         f = re.search('_\d+.\d+[pm]*_',html)
#         ri = f.group()[1:-1]
#         # Reading the file 
#         text = html_file.read() 
#         rel_section = re.search('Event[^{]+'+str(junction)+'[\s\S]+}',text)
#         while not rel_section:
#             junction = j_to_j[junction]
#             rel_section = re.search('Event[^{]+'+str(junction)+'[\s\S]+}',text)
#         response = rel_section.group()
#         evnt = re.search('(?<=Event\": \")[^\"]+',rel_section.group()) 
#         #print(huc,junction,evnt.group())
#         events_dict[huc][original_j][ri] = evnt.group()
#     return events_dict

def get_sst_storms_by_recurrence_us_huc(events_dict,huc,htmls,junction):
    #redefined get_sst_storms_by_recurrence_us_huc to not grab any events in the next junction
    original_j = junction
    for html in htmls:
        html_file = open(html, "r") 
        f = re.search('_\d+\.\d+[pm]*_',html)
        ri = f.group()[1:-1]
        # Reading the file 
        text = html_file.read() 
        rel_section = re.search('Event[^{]+'+str(junction)+'[\s\S]+}',text)
        
        # #special case to find downstream junction of us huc if not 
        # if not rel_section:
        #     #print(f'{junction} not in {huc} html file. Finding most downstream junction in huc10')
        #     ds_junction = j_to_j[junction]
        #     junction_huc = re.search('HUC_[0-9]{3}',junction).group()
        #     ds_junction_huc = re.search('HUC_[0-9]{3}',ds_junction).group()
        #     assert junction_huc != ds_junction_huc, f'downstream junction/sink/res of {junction} does not appear to actually be the most downstream'
        #     for key, val in j_to_j.items():
        #         if key.find(junction_huc+'_J') >=0 and val.find(ds_junction_huc+'_J') >=0:
        #             junction = key
        #             #print(f'new ds junction is {key}')
        #         else:
        #             pass
        #     rel_section = re.search('Event[^{]+'+str(junction)+'[\s\S]+}',text)
        assert rel_section, f"logic in code is not foolproof. Talk to Curtis: {rel_section}"
        response = rel_section.group()
        evnt = re.search('(?<=Event\": \")[^\"]+',rel_section.group()) 
        #print(huc,junction,evnt.group())
        events_dict[huc][original_j][ri] = evnt.group()
        
    return events_dict

def update_oos_bcs(in_geom_path,index,bc_dict,domain_geo, geo_out):
    #read old text
    with open(in_geom_path, "r+") as f:
        file_contents = f.read()
    new_text = ''
    #loop through so nothing happens if no updates are needed
    for i in index:
        #add bc lines (which are offset)
        geom = bc_dict['geometry'][i]
        first_coord = ' , '.join(str(coord) for coord in geom.coords[0])
        mid_coord = ' , '.join(str(coord) for coord in geom.line_interpolate_point(geom.length/2).coords[0])
        last_coord = ' , '.join(str(coord) for coord in geom.coords[-1])
        coord_len = len(geom.coords)
        junc = bc_dict['junction'][i]
        ext_bc_template = [f'BC Line Name={junc}',
                f'BC Line Storage Area= {domain_geo}',
                f'BC Line Start Position={first_coord}', 
                f'BC Line Middle Position={mid_coord}', 
                f'BC Line End Position={last_coord} ',
                f'BC Line Arc= {coord_len}']
        for j, coord in enumerate(geom.coords):
            if j % 2 == 0:
                coord_str = f"{coord[0]:16f}"+f"{coord[1]:16f}"
            else:
                coord_str += f"{coord[0]:16f}"+f"{coord[1]:16f}"
                ext_bc_template.append(coord_str)
                coord_str = ''
        if coord_str != '':
            ext_bc_template.append(coord_str)
        ext_bc_template.append('BC Line Text Position= 1.79769313486232E+308 , 1.79769313486232E+308')
        new_text += '\n'.join(ext_bc_template)
        new_text += '\n'

        #find insertion points
        bc_name_exist = bc_dict['huc'][i]
        geo__oos_bc_finder= f'BC Line Name=inflow_from_{bc_name_exist}[\d\s\S]+?(?=BC Line Name)|(?=LCMann Time)'
        bc_insertion = re.search(geo__oos_bc_finder,file_contents)
        start = file_contents[:bc_insertion.start()]
        middle = new_text
        end = file_contents[bc_insertion.end():]
        geo_out_str = start+middle+end

        with open(geo_out, "w") as g_o:
            g_o.write(geo_out_str)


#######################################################################################
#combine the schematics for all GDB in the project area. potentially do this once? 
def locate_string_gdb_to_concat_gdf(gdb_path_list,locater_string,coordinate_sys=None):
    """Creates a geodataframe by searching through a list of geodatabases/geopackages for a specific layer name.

    Args:
        gdb_path_list (list): list of paths to geodatabases/geopackages
        locater_string (string): string to identify layer name within geodatabases/geopackages
        coordinate_sys (string, optional): EPSG, e.g. "4326". Defaults to None.

    Returns:
        Geodataframe: _description_
    """
    gdb_quantity = len(gdb_path_list)
    
    if gdb_quantity == 0:
        print(f'cannot identify any {locater_string} from the path list provided')
        gdf_x = gpd.GeoDataFrame()
        gdf_x['geometry'] = None
        coord = coordinate_sys
    
    if gdb_quantity != 0:
        layers_list = fiona.listlayers(gdb_path_list[0])
        if locater_string in layers_list:
            gdf_x = gpd.read_file(gdb_path_list[0], layer=f"{locater_string}")
            coord = gdf_x.crs
            gdf_x['Source_Path'] = os.path.basename(gdb_path_list[0])
            if coordinate_sys!= None:
                coord = coordinate_sys
                gdf_x = gdf_x.to_crs(coord)
            else:
                pass
        else:
            gdf_x = gpd.GeoDataFrame()
            gdf_x['geometry'] = None
            coord = None
    
    if gdb_quantity > 1:
        for i in range(1,gdb_quantity):
            layers_list2 = fiona.listlayers(gdb_path_list[i])
            if locater_string in layers_list2:
                gdf_x2 = gpd.read_file(gdb_path_list[i], layer=f"{locater_string}")
                coord2 = gdf_x2.crs
                gdf_x2['Source_Path'] = os.path.basename(gdb_path_list[i])
                if coord == None:
                    coord = coord2
                if coord2 != coord:
                    print("\ncoordinate systems are not the same: will convert to match first gdb/gpkg in path list")
                    gdf_x2 = gdf_x2.to_crs(gdf_x.crs)
                gdf_x = pd.concat([gdf_x,gdf_x2])
                gdf_x.reset_index(drop=True, inplace=True)
            else:
                pass
    
    print(f'After searching through {gdb_quantity} geodatabases/geopackages, there are now {len(gdf_x)} \"{locater_string}\" features')
    return gdf_x


