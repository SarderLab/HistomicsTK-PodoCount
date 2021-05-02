"""
Created on Mon Nov 23 9:00:07 2020
@author: Briana Santo

"""
import numpy as np
import skimage as sk
import scipy as sp
import pandas as pd
import matplotlib.pyplot as mplp
import os
import cv2
import openslide
import glob
import time
import warnings
from xml_to_mask import xml_to_mask
from get_wsi_mask import get_wsi_mask
from get_boundary import get_boundary
from get_glom_props import get_glom_props
from get_pod_feat_spaces import get_pod_feat_spaces
from stain_decon import stain_decon
from get_pod_props import get_pod_props
from skimage import color, morphology
from skimage.transform import resize
import argparse

#filter warnings
warnings.filterwarnings("ignore")

#get current working directory
'''
python3 ../BAS_folder/pod_quantification_main.py -A '{}' -B '{}' -D {} -P '{}' -L '{}' -M '{}'".format(args2.inputImageFile, args2.inputAnnotationFile1,args2.inputAnnotationFile2,args2.slider, args2.ihc_gauss_sd, args2.num_sections, args2.outputAnnotationFile1 ,args2.outputAnnotationFile2 , args2.csvFile)
python3 pod_quantification_main.py -A '{}' -B '{}' -D {} -P '{}' -L '{}' -M '{}'".format(args2.inputImageFile, args2.inputAnnotationFile1,args2.inputAnnotationFile2,args2.slider, args2.ihc_gauss_sd, args2.num_sections, args2.outputAnnotationFile1 ,args2.outputAnnotationFile2 , args2.csvFile)
python3 pod_quant_main_serv.py -A '1.svs' -B '1glom.xml' -D 2.5 -P 'out1.xml' -L 'out2.xml' -M 'out3.csv'

'''
parser = argparse.ArgumentParser(description = '')
parser.add_argument('-A','--inputsvs',type = str, metavar = '',required = True,help = 'image name')
parser.add_argument('-B','--glomxml',type = str, metavar = '',required = True,help = ' glom xml')
parser.add_argument('-D','--var1',type = float, metavar = '',required = True,help = 'slider')
parser.add_argument('-P','--outxml1',type = str, metavar = '',required = True,help = 'outxml1')
parser.add_argument('-L','--outxml2',type = str, metavar = '',required = True,help = 'outxml2')
parser.add_argument('-M','--outcsv',type = str, metavar = '',required = True,help = 'outcsv')

args = parser.parse_args()


WSI_file = args.inputsvs
WSI_glom_xml = args.glomxml
slider = args.var1
xml_contour = args.outxml1
xml_counter = args.outxml2
csv_path = args.outcsv

#WSI general info
#ftype = '.svs'
#WSIs = str(cwd + '/WSIs/*' + ftype)
#glom_xmls = str(cwd + '/glom_xmls/*.xml')

#WSI_dir = glob.glob(WSIs)
#glom_xmls_dir = glob.glob(glom_xmls)

# if ftype is .ndpi then

#Parameters
#slider = 2.25
num_sections = 1

#Main Script
for WSI in range(1,2):
    start_time = time.time()

#    WSI_file = WSI_dir[WSI]
    ftype = WSI_file[-4:]
    WSI_name = WSI_file.split('/')
    WSI_name = WSI_name[-1]
    WSI_name = WSI_name.split('.')
    WSI_name = WSI_name[0]
    print('\n')
    print('--- Working on: '+str(WSI_name)+' ---\n')
    WSI_file = openslide.open_slide(WSI_file)
    WSI_meta = (float(WSI_file.properties[openslide.PROPERTY_NAME_MPP_X])+float(WSI_file.properties[openslide.PROPERTY_NAME_MPP_Y]))/2
    dist_mpp, area_mpp2 = WSI_meta, WSI_meta**2
    WSI_cols,WSI_rows = WSI_file.dimensions[0],WSI_file.dimensions[1]
    
#    ftype = WSI_file[-4:]
    
    if ftype == '.ndpi':
        level = 2
    else:
        level = 1
    WSI_levels = WSI_file.level_dimensions[level]
    downsample_factor=16
    df2 = np.sqrt(downsample_factor)

    #get whole-slide glom mask
#    WSI_glom_xml = cwd + '/glom_xmls/' + WSI_name + '.xml'
    WSI_glom_mask = xml_to_mask(WSI_glom_xml, (0,0), (WSI_cols,WSI_rows), downsample_factor=downsample_factor, verbose=0)
    WSI_glom_mask = np.array(WSI_glom_mask)

    #get whole-slide mask
    WSI_downsample = np.array(WSI_file.read_region((0,0),level,(WSI_levels[0],WSI_levels[1])),dtype = 'uint8')
    WSI_downsample = get_wsi_mask(WSI_downsample[:,:,0:3],WSI_glom_mask,num_sections)

    #xml files - initiation
    #xml_counter = WSI_name +'_counter.xml'
    #xml_contour = WSI_name + '_contour.xml'

    xml_counter = open(xml_counter,'w')
    xml_contour = open(xml_contour,'w')

    xml_counter.write('<Annotations>\n\t<Annotation Id="1">\n\t\t<Attributes>\n\t\t\t<Attribute Id="0" Name="Pod" Value="" />\n\t\t</Attributes>\n\t\t<Regions>')
    xml_contour.write('<Annotations>\n\t<Annotation Id="1">\n\t\t<Attributes>\n\t\t\t<Attribute Id="0" Name="Pod" Value="" />\n\t\t</Attributes>\n\t\t<Regions>')

    count = 0

    #get ROI coordinates
    print('-- Step 1: Glomerulus localization and quantification --\n')
    bbs, total_gloms, glom_feat_labels, glom_feat_qty, glom_feat_array = get_glom_props(WSI_glom_mask,WSI_downsample,num_sections,dist_mpp,area_mpp2,df2)

    #define pod feature spaces
    pod_feat_labels, pod_feat_qty, pod_feat_array = get_pod_feat_spaces(total_gloms)

    print('-- Step 2: Podocyte detection and quantification --\n')
    for bb in range(len(bbs)):
        bb_iter = bb
        bb = bbs[bb]
        x_start,y_start,x_stop,y_stop = int(df2*bb[0]),int(df2*bb[1]),int(df2*bb[2]),int(df2*bb[3])
        x_length = x_stop-x_start
        y_length = y_stop-y_start

        roi = np.array(WSI_file.read_region((y_start,x_start),0,(y_length,x_length)),dtype = 'uint8')
        roi = roi[:,:,0:3]
        rows,cols,dims = roi.shape

        glom_mask = WSI_glom_mask[bb[0]:bb[2],bb[1]:bb[3]]
    #    glom_mask = sp.misc.imresize(glom_mask,[rows,cols],interp='bilinear', mode=None)
        glom_mask = (resize(glom_mask,[rows,cols],anti_aliasing=True)>0.01)*1

        xml_counter, xml_contour, count, pod_feat_vector = get_pod_props(roi,glom_mask,slider,x_start,y_start,xml_counter,xml_contour,count,dist_mpp,area_mpp2)
        pod_feat_array[:,bb_iter] = pod_feat_vector

    #xml files - finalization
    xml_counter.write('\n\t\t</Regions>\n\t</Annotation>\n</Annotations>')
    xml_contour.write('\n\t\t</Regions>\n\t</Annotation>\n</Annotations>')

    xml_counter.close()
    xml_contour.close()

    print('-- Step 3: Feature file creation --\n')
    final_feat_array = np.vstack([glom_feat_array.T,pod_feat_array])
    final_feat_labels = glom_feat_labels + pod_feat_labels
    final_feat_DF = pd.DataFrame(final_feat_array,index=final_feat_labels)
    #csv_path = WSI_name + '_Features.csv'
    final_feat_DF.to_csv(csv_path,index=True,columns=None)

    print('--- Completed: '+str(WSI_name)+' ---\n')
    end_time = time.time() - start_time
    print("--- %s seconds for whole-slide analysis ---" % (end_time))

print('\n')
print('--- Completed full cohort ---')
print('\n')
