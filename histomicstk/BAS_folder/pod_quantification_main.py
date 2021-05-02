"""
Created on Fri Sept 4 11:40:33 2020
@author: Briana Santo
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import skimage as sk
import scipy as sp
import pandas as pd

import sys
sys.path.append("..")

import os.path
#import matplotlib.pyplot as plt
import cv2
#import glob
#import os
import openslide
from BAS_folder import xml_to_mask
from BAS_folder import get_glom_props
from BAS_folder import get_pod_feat_spaces
from BAS_folder import get_pod_props
from skimage import color,morphology
from skimage.transform import resize
#from skimage import segmentation
#from skimage.measure import find_contours
import argparse


parser = argparse.ArgumentParser(description = '')
parser.add_argument('-A','--inputsvs',type = str, metavar = '',required = True,help = 'image name')
parser.add_argument('-B','--glomxml',type = str, metavar = '',required = True,help = ' glom xml')
parser.add_argument('-C','--medxml',type = str, metavar = '',required = True,help = 'med xml')
parser.add_argument('-D','--var1',type = float, metavar = '',required = True,help = 'slider')
parser.add_argument('-E','--var2',type = int, metavar = '',required = True,help = 'ihc_gauss_sd')
parser.add_argument('-F','--var3',type = int, metavar = '',required = True,help = 'num_sections')
parser.add_argument('-P','--outxml1',type = str, metavar = '',required = True,help = 'outxml1')
parser.add_argument('-L','--outxml2',type = str, metavar = '',required = True,help = 'outxml2')
parser.add_argument('-M','--outcsv',type = str, metavar = '',required = True,help = 'outcsv')

args = parser.parse_args()


WSI_file = args.inputsvs
WSI_glom_xml = args.glomxml
WSI_medulla_xml = args.medxml
slider = args.var1
ihc_gauss_sd = args.var2
num_sections = args.var3
xml_contour = args.outxml1
xml_counter = args.outxml2
csv_path = args.outcsv

#ContourXML = WSIfile, split, + '.xml'


#NOTE TO DG ON INPUTS AND OUTPUTS#
#INPUTS are: WSI, GLOM XML, MEDULLA XML
#OUTPUTS are: Contour XML (named 'WSI_name' + '.xml'), Counter XML, FEATURE CSV FILE
#PARAMETERS THAT NEED TO BE INPUT ARE:

#Parameters
#slider = 2.5
#ihc_gauss_sd = 2
#num_sections = 1

#Functions
def get_boundary(image):
    boundary_pts,_ = cv2.findContours(np.uint8(image), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    boundary_pts = np.vstack(boundary_pts)
    boundary_pts = boundary_pts.reshape(len(boundary_pts),2)
    return boundary_pts

def get_wsi_mask(image):
    wsi_mask = color.rgb2gray(image)
    wsi_mask = wsi_mask<(wsi_mask.mean())
    wsi_mask = sp.ndimage.morphology.binary_fill_holes(wsi_mask)
    wsi_labels, num_labels = sp.ndimage.label(wsi_mask)
    wsi_props = sk.measure.regionprops(wsi_labels)
    areas_labels = np.empty([num_labels,2])
    for label in range(num_labels):
        areas_labels[label,0] = wsi_props[label].area
        areas_labels[label,1] = (label+1)
    areas_labels = areas_labels[areas_labels[:,0].argsort()]
    wsi_mask = np.zeros(wsi_mask.shape)
    for section in range(num_sections):
        label = areas_labels[-(section+1),1]
        wsi_mask[wsi_labels==label] = 1
    return(wsi_mask)

#Main Script
#WSI_file = [INPUT WSI NAME]
WSI_name = WSI_file.split('.')
WSI_name = WSI_name[0]
WSI_ext = WSI_name[1]
print('Working on:',WSI_name)
WSI_file = openslide.open_slide(WSI_file)
WSI_meta = (float(WSI_file.properties[openslide.PROPERTY_NAME_MPP_X])+float(WSI_file.properties[openslide.PROPERTY_NAME_MPP_Y]))/2
dist_mpp, area_mpp2 = WSI_meta, WSI_meta**2
WSI_cols,WSI_rows = WSI_file.dimensions[0],WSI_file.dimensions[1]
WSI_levels = WSI_file.level_dimensions[1]
downsample_factor=16
df2 = np.sqrt(downsample_factor)

#get whole-slide mask
WSI_downsample = np.array(WSI_file.read_region((0,0),1,(WSI_levels[0],WSI_levels[1])),dtype = 'uint8')
WSI_downsample = get_wsi_mask(WSI_downsample[:,:,0:3])

#get whole-slide glom mask
#WSI_glom_xml = [INPUT GLOM XML NAME]
WSI_glom_mask = xml_to_mask.xml_to_mask(WSI_glom_xml, (0,0), (WSI_cols,WSI_rows), downsample_factor=downsample_factor, verbose=0)
WSI_glom_mask = np.array(WSI_glom_mask)

#get whole-slide medulla mask
#WSI_medulla_xml = [INPUT MEDULLA XML NAME]
WSI_medulla_mask = xml_to_mask.xml_to_mask(WSI_medulla_xml, (0,0), (WSI_cols,WSI_rows), downsample_factor=16, verbose=0)
WSI_medulla_mask = np.array(WSI_medulla_mask)

#xml files - initiation
#xml_counter = 'counter.xml'
#xml_contour = WSI_name + '.xml'

xml_counter = open(xml_counter,'w')
xml_contour = open(xml_contour,'w')

xml_counter.write('<Annotations>\n\t<Annotation Id="1">\n\t\t<Attributes>\n\t\t\t<Attribute Id="0" Name="Glom" Value="" />\n\t\t</Attributes>\n\t\t<Regions>')
xml_contour.write('<Annotations>\n\t<Annotation Id="1">\n\t\t<Attributes>\n\t\t\t<Attribute Id="0" Name="Glom" Value="" />\n\t\t</Attributes>\n\t\t<Regions>')

count = 0

#get ROI coordinates
print('Step 1: Glomerulus localization and quantification')
bbs, total_gloms, glom_feat_labels, glom_feat_qty, glom_feat_array = get_glom_props.get_glom_props(WSI_glom_mask,WSI_downsample,WSI_medulla_mask,num_sections,dist_mpp,area_mpp2,df2)

#define pod feature spaces
pod_feat_labels, pod_feat_qty, pod_feat_array = get_pod_feat_spaces.get_pod_feat_spaces(total_gloms)

print('Step 2: Podocyte detection and quantification')
for bb in range(len(bbs)):
    bb_iter = bb
    bb = bbs[bb]
    x_start,y_start,x_stop,y_stop = int(df2*bb[0]),int(df2*bb[1]),int(df2*bb[2]),int(df2*bb[3])
    x_length = x_stop-x_start
    y_length = y_stop-y_start

    ROI = np.array(WSI_file.read_region((y_start,x_start),0,(y_length,x_length)),dtype = 'uint8')
    ROI = ROI[:,:,0:3]
    rows,cols,dims = ROI.shape

    glom_mask = WSI_glom_mask[bb[0]:bb[2],bb[1]:bb[3]]
#    glom_mask = sp.misc.imresize(glom_mask,[rows,cols],interp='bilinear', mode=None)
    glom_mask = (resize(glom_mask,[rows,cols],anti_aliasing=True)>0.01)*1
    
#    plt.figure()
#    plt.imshow(glom_mask)
#    plt.show()
    
    xml_counter, xml_contour, count, pod_feat_vector = get_pod_props.get_pod_props(ROI,glom_mask,slider,ihc_gauss_sd,x_start,y_start,xml_counter,xml_contour,count,dist_mpp,area_mpp2)
    pod_feat_array[:,bb_iter] = pod_feat_vector

#xml files - finalization
xml_counter.write('\n\t\t</Regions>\n\t</Annotation>\n</Annotations>')
xml_contour.write('\n\t\t</Regions>\n\t</Annotation>\n</Annotations>')

xml_counter.close()
xml_contour.close()

print('Step 3: Feature file creation')
final_feat_array = np.vstack([glom_feat_array.T,pod_feat_array])
final_feat_labels = glom_feat_labels + pod_feat_labels
final_feat_DF = pd.DataFrame(final_feat_array,index=final_feat_labels)
#csv_path = WSI_name + '_Features.csv'
final_feat_DF.to_csv(csv_path,index=True,columns=None)

print('Completed:',WSI_name)
