"""
Created on Fri Sept 4 11:40:33 2020
@author: Briana Santo
"""
import numpy as np
import skimage as sk
import scipy as sp
#import matplotlib.pyplot as mplp
import cv2
#from skimage import color,morphology
#from skimage import segmentation
#from skimage.measure import find_contours

def get_boundary(image):
    boundary_pts,_ = cv2.findContours(np.uint8(image), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    boundary_pts = np.vstack(boundary_pts)
    boundary_pts = boundary_pts.reshape(len(boundary_pts),2)
    return boundary_pts

def get_zones(glom_count,centroids,cortex_boundary,medulla_boundary,dist_mpp,area_mpp2,df2):
    cortex_dists = []
    medulla_dists = []
    for glom in range(glom_count):
        centroid = np.array(centroids[glom]).reshape(1,-1)
        distances = sp.spatial.distance.cdist(centroid, cortex_boundary, 'euclidean')
        cortex_dists.append((np.min(distances))*df2*dist_mpp)

        centroid = np.array(centroids[glom]).reshape(1,-1)
        distances = sp.spatial.distance.cdist(centroid, medulla_boundary, 'euclidean')
        medulla_dists.append((np.min(distances))*df2*dist_mpp)

    cortex_dists = np.array(cortex_dists)
    medulla_dists = np.array(medulla_dists)
    depths = medulla_dists/(np.add(cortex_dists,medulla_dists))

    _3zoneID = np.empty([len(depths),1])
    _5zoneID = np.empty([len(depths),1])
    for depth in range(len(depths)):
        if depths[depth]<=0.4:
            _3zoneID[depth] = 1
        elif 0.4<depths[depth]<=0.7:
            _3zoneID[depth] = 2
        else: _3zoneID[depth] = 3

    for depth in range(len(depths)):
        if depths[depth]<=0.4:
            _5zoneID[depth] = 1
        elif 0.4<depths[depth]<=0.55:
            _5zoneID[depth] = 2
        elif 0.55<depths[depth]<=0.7:
            _5zoneID[depth] = 3
        elif 0.7<depths[depth]<=0.85:
            _5zoneID[depth] = 4
        else: _5zoneID[depth] = 5

    return(cortex_dists,medulla_dists,_3zoneID,_5zoneID)

def get_spatial_density(glom_count,centroids,dist_mpp,area_mpp2,df2):
    inter_glom_dists = []
    for glom in range(glom_count):
        centroid = np.array(centroids[glom]).reshape(1,-1)
        distances = sp.spatial.distance.cdist(centroid, centroids, 'euclidean')
        distances = np.sort(distances)
        inter_glom_dists.append((np.sum(distances[0:5]))*df2*dist_mpp)

    glom_density = (1.0)/np.array(inter_glom_dists)
    return(inter_glom_dists,glom_density)

###GLOM FEATURE ENGINEERING AND EXTRACTION###
def get_glom_props(glom_image,tissue_image,medulla_image,num_sections,dist_mpp,area_mpp2,df2):
    if num_sections >1:
        sections_label, section_count= sp.ndimage.label(tissue_image)
        s1_gloms = np.logical_and(sections_label==1,glom_image)
        s2_gloms = np.logical_and(sections_label==2,glom_image)

        s1_glom_labels,s1_glom_count = sp.ndimage.label(s1_gloms)
        s2_glom_labels,s2_glom_count = sp.ndimage.label(s2_gloms)
        total_gloms = s1_glom_count + s2_glom_count
        s1_gen_props = sk.measure.regionprops(s1_glom_labels)

        glom_feat_array = []
        s_label = []
        glomID = []
        areas = []
        bbs = []
        bb_areas = []
        centroids = []
        convex_areas = []
        eccentricities = []
        equiv_diams = []
        extents = []
        major_axis_lengths = []
        minor_axis_lengths = []
        orientations = []
        perimeters = []
        solidities = []

        for glom in range(s1_glom_count):

            s_label.append(1)
            glomID.append(glom+1)
            areas.append((s1_gen_props[glom].area)*df2*area_mpp2)
            bbs.append(s1_gen_props[glom].bbox)
            bb_areas.append((s1_gen_props[glom].bbox_area)*df2*area_mpp2)
            centroids.append(s1_gen_props[glom].centroid)
            convex_areas.append((s1_gen_props[glom].convex_area)*df2*area_mpp2)
            eccentricities.append(s1_gen_props[glom].eccentricity)
            equiv_diams.append((s1_gen_props[glom].equivalent_diameter)*df2*dist_mpp)
            extents.append(s1_gen_props[glom].extent)
            major_axis_lengths.append((s1_gen_props[glom].major_axis_length)*df2*dist_mpp)
            minor_axis_lengths.append((s1_gen_props[glom].minor_axis_length)*df2*dist_mpp)
            orientations.append(s1_gen_props[glom].orientation)
            perimeters.append((s1_gen_props[glom].perimeter)*df2*dist_mpp)
            solidities.append(s1_gen_props[glom].solidity)

        #find glom spatial density
        [s1_inter_glom_dists,s1_glom_density] = get_spatial_density(s1_glom_count,centroids,dist_mpp,area_mpp2,df2)

        #find superficial cortical and medullary distances
        s1_med = np.logical_and(sections_label==1,medulla_image)
        s2_med = np.logical_and(sections_label==2,medulla_image)

        s1_cortex_boundary = get_boundary(sections_label==1)
        s1_medulla_boundary = get_boundary(s1_med)
        s2_cortex_boundary = get_boundary(sections_label==2)
        s2_medulla_boundary = get_boundary(s2_med)

        s1_cortex_boundary[:, [1, 0]] = s1_cortex_boundary[:, [0, 1]]
        s1_medulla_boundary[:, [1, 0]] = s1_medulla_boundary[:, [0, 1]]
        s2_cortex_boundary[:, [1, 0]] = s2_cortex_boundary[:, [0, 1]]
        s2_medulla_boundary[:, [1, 0]] = s2_medulla_boundary[:, [0, 1]]

        [s1_cortex_dists,s1_medulla_dists,s1_3zoneID,s1_5zoneID] = get_zones(s1_glom_count,centroids,s1_cortex_boundary,s1_medulla_boundary,dist_mpp,area_mpp2,df2)
        #print('stack:',np.hstack([np.array(glomID).reshape([73,1]),np.array(areas).reshape([73,1]),s1_3zoneID,s1_5zoneID]))
        #mplp.imshow(s1_glom_labels)
        #mplp.show()

        s2_gen_props = sk.measure.regionprops(s2_glom_labels)
        for glom in range(s2_glom_count):
            s_label.append(2)
            glomID.append(glom+1)
            areas.append((s2_gen_props[glom].area)*df2*area_mpp2)
            bbs.append(s2_gen_props[glom].bbox)
            bb_areas.append((s2_gen_props[glom].bbox_area)*df2*area_mpp2)
            centroids.append(s2_gen_props[glom].centroid)
            convex_areas.append((s2_gen_props[glom].convex_area)*df2*area_mpp2)
            eccentricities.append(s2_gen_props[glom].eccentricity)
            equiv_diams.append((s2_gen_props[glom].equivalent_diameter)*df2*dist_mpp)
            extents.append(s2_gen_props[glom].extent)
            major_axis_lengths.append((s2_gen_props[glom].major_axis_length)*df2*dist_mpp)
            minor_axis_lengths.append((s2_gen_props[glom].minor_axis_length)*df2*dist_mpp)
            orientations.append(s2_gen_props[glom].orientation)
            perimeters.append((s2_gen_props[glom].perimeter)*df2*dist_mpp)
            solidities.append(s2_gen_props[glom].solidity)

        [s2_inter_glom_dists,s2_glom_density] = get_spatial_density(s2_glom_count,centroids,dist_mpp,area_mpp2,df2)
        [s2_cortex_dists,s2_medulla_dists,s2_3zoneID,s2_5zoneID] = get_zones(s2_glom_count,centroids,s2_cortex_boundary,s2_medulla_boundary,dist_mpp,area_mpp2,df2)

        inter_glom_dists = np.vstack([np.array(s1_inter_glom_dists).reshape(s1_glom_count,1),np.array(s2_inter_glom_dists).reshape(s2_glom_count,1)])
        glom_densities = np.vstack([np.array(s1_glom_density).reshape(s1_glom_count,1),np.array(s2_glom_density).reshape(s2_glom_count,1)])
        cortex_dists = np.vstack([np.array(s1_cortex_dists).reshape(s1_glom_count,1),np.array(s2_cortex_dists).reshape(s2_glom_count,1)])
        medulla_dists = np.vstack([np.array(s1_medulla_dists).reshape(s1_glom_count,1),np.array(s2_medulla_dists).reshape(s2_glom_count,1)])
        _3zoneIDs = np.vstack([np.array(s1_3zoneID).reshape(s1_glom_count,1),np.array(s2_3zoneID).reshape(s2_glom_count,1)])
        _5zoneIDs = np.vstack([np.array(s1_5zoneID).reshape(s1_glom_count,1),np.array(s2_5zoneID).reshape(s2_glom_count,1)])

    elif num_sections == 1:
        glom_label, glom_count = sp.ndimage.label(glom_image)
        gen_props = sk.measure.regionprops(glom_label)
        total_gloms = glom_count

        glom_feat_array = []
        s_label = []
        glomID = []
        areas = []
        bbs = []
        bb_areas = []
        centroids = []
        convex_areas = []
        eccentricities = []
        equiv_diams = []
        extents = []
        major_axis_lengths = []
        minor_axis_lengths = []
        orientations = []
        perimeters = []
        solidities = []

        for glom in range(glom_count):
            s_label.append(1)
            glomID.append(glom+1)
            areas.append((gen_props[glom].area)*df2*area_mpp2)
            bbs.append(gen_props[glom].bbox)
            bb_areas.append((gen_props[glom].bbox_area)*df2*area_mpp2)
            centroids.append(gen_props[glom].centroid)
            convex_areas.append((gen_props[glom].convex_area)*df2*area_mpp2)
            eccentricities.append(gen_props[glom].eccentricity)
            equiv_diams.append((gen_props[glom].equivalent_diameter)*df2*dist_mpp)
            extents.append(gen_props[glom].extent)
            major_axis_lengths.append((gen_props[glom].major_axis_length)*df2*dist_mpp)
            minor_axis_lengths.append((gen_props[glom].minor_axis_length)*df2*dist_mpp)
            orientations.append(gen_props[glom].orientation)
            perimeters.append((gen_props[glom].perimeter)*df2*dist_mpp)
            solidities.append(gen_props[glom].solidity)

        #find glom spatial density
        [inter_glom_dists,glom_densities] = get_spatial_density(glom_count,centroids,dist_mpp,area_mpp2,df2)

        #find superficial cortical and medullary distances
        cortex_boundary = get_boundary(tissue_image)
        medulla_boundary = get_boundary(medulla_image)
        cortex_boundary[:, [1, 0]] = cortex_boundary[:, [0, 1]]
        medulla_boundary[:, [1, 0]] = medulla_boundary[:, [0, 1]]

        [cortex_dists,medulla_dists,_3zoneIDs,_5zoneIDs] = get_zones(glom_count,centroids,cortex_boundary,medulla_boundary,dist_mpp,area_mpp2,df2)

    glom_feat_labels = ['section_label','glomID','glom_areas','glom_bb_areas','glom_convex_areas','glom_eccentricities','glom_equiv_diams','glom_extents','glom_major_axis_lengths','glom_minor_axis_lengths','glom_orientations','glom_perimeters','glom_solidities','inter_glom_dists','glom_density','glom_cortex_dists','glom_medulla_dists','glom_3zoneID','glom_5zoneID']
    glom_feat_qty = len(glom_feat_labels)

    glom_feat_array =( np.hstack([np.array(s_label).reshape([total_gloms,1]),np.array(glomID).reshape([total_gloms,1]),
    np.array(areas).reshape([total_gloms,1]),np.array(bb_areas).reshape([total_gloms,1]),np.array(convex_areas).reshape([total_gloms,1]),
    np.array(eccentricities).reshape([total_gloms,1]),np.array(equiv_diams).reshape([total_gloms,1]),np.array(extents).reshape([total_gloms,1]),
    np.array(major_axis_lengths).reshape([total_gloms,1]),np.array(minor_axis_lengths).reshape([total_gloms,1]),np.array(orientations).reshape([total_gloms,1]),
    np.array(perimeters).reshape([total_gloms,1]),np.array(solidities).reshape([total_gloms,1]),np.array(inter_glom_dists).reshape([total_gloms,1]),
    np.array(glom_densities).reshape([total_gloms,1]),np.array(cortex_dists).reshape([total_gloms,1]),np.array(medulla_dists).reshape([total_gloms,1]),np.array(_3zoneIDs).reshape([total_gloms,1]),np.array(_5zoneIDs).reshape([total_gloms,1])]) )

    return np.array(bbs), total_gloms, glom_feat_labels, glom_feat_qty, glom_feat_array


    print('features done')
