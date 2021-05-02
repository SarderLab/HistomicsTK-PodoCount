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
from skimage import segmentation
#from skimage.measure import find_contours

#Functions
def stain_decon(image):
    stain1 = np.array([0.3697298, 0.61498046, 0.69649047])
    stain2 = np.array([0.4091243, 0.8440652, 0.34665722])
    Residual = np.array([0.47995895, 0.6926196, 0.5384398])
    Res = np.cross(stain1, stain2)
    hdab_rgb = np.array([stain1/np.linalg.norm(stain1), stain2/np.linalg.norm(stain2), Res/np.linalg.norm(Res)])
    hdab_rgb = hdab_rgb.T
    rgb_hdab = np.linalg.inv(hdab_rgb);
    sep_stains = sk.color.separate_stains(image,rgb_hdab)
    return sep_stains

def get_boundary(image):
    boundary_pts,_ = cv2.findContours(np.uint8(image), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    boundary_pts = np.vstack(boundary_pts)
    boundary_pts = boundary_pts.reshape(len(boundary_pts),2)
    return boundary_pts

def get_pod_props(ROI,glom_mask,slider,ihc_gauss_sd,x_start,y_start,xml_counter,xml_contour,count,dist_mpp,area_mpp2):
    #Parameters
    dt_gauss_sd = 1
    min_area = 150
    emt_thresh = 0.1

    ###PODOCYTE SEGMENTATION###
    podo_decon = stain_decon(ROI)
    ihc = podo_decon[:,:,2]
    lumen = podo_decon[:,:,1]
    pas = podo_decon[:,:,0]

    ihc = np.uint8((sk.filters.gaussian(ihc,ihc_gauss_sd)*255))
    ihc_mean = ihc.mean()
    ihc_std = ihc.std()
    ihc_thresh = (ihc_mean + (slider*ihc_std))

    nuclei = (ihc>ihc_thresh)
    nuclei = sp.ndimage.morphology.binary_fill_holes(nuclei)
    nuclei_uni = sk.morphology.remove_small_objects(nuclei, min_size=100, connectivity=8)

    nuclei_dt = sp.ndimage.morphology.distance_transform_edt(~nuclei)
    nuclei_dt = sk.filters.gaussian(nuclei_dt,dt_gauss_sd)
    nuclei_dt_max = sk.morphology.h_minima(nuclei_dt,emt_thresh*nuclei_dt.max())

    nuclei_label, num_labels = sp.ndimage.label(nuclei_uni - nuclei_dt_max)
    nuclei_props = sk.measure.regionprops(nuclei_label)
    nuc_euls = []

    for nuc in range(num_labels):
        nuc_eul = nuclei_props[nuc].euler_number
        if nuc_eul <0:
            nuc_euls.append(nuc)

    nuclei_multi = np.zeros(nuclei_label.shape)
    if np.sum(nuc_euls)>0:
        for nuc in range(len(nuc_euls)):
            label = nuc_euls[nuc] + 1
            nuclei_multi[nuclei_label==label] = 1
            nuclei_uni[nuclei_label==label] = 0

        nuclei_multi = sp.ndimage.morphology.binary_fill_holes(nuclei_multi)
        nuclei_multi_dt = sp.ndimage.morphology.distance_transform_edt(nuclei_multi)
        markers = np.logical_and(nuclei_dt_max,nuclei_multi)
        nuclei_multi_wt = segmentation.watershed(-nuclei_multi_dt,markers,watershed_line=True)

        separated_podocytes = np.logical_and(nuclei_multi,(nuclei_multi_wt==0))>0
        separated_podocytes = sk.morphology.remove_small_objects(separated_podocytes, min_size=100, connectivity=8)
        separated_podocytes = np.logical_xor(separated_podocytes,nuclei_uni)

    else:
        separated_podocytes = nuclei_uni
        markers = np.zeros(nuclei_label.shape)

    rows,cols = np.shape(separated_podocytes)
    podocytes = separated_podocytes*(np.resize(glom_mask,(rows,cols)))

    ###PODOCYTE FEATURE ENGINEERING AND EXTRACTION###
    podocyte_label, podocyte_count = sp.ndimage.label(podocytes)

    if podocyte_count>0:
        gen_props = sk.measure.regionprops(podocyte_label,ihc)

        pod_feat_vector = []
        areas = []
        centroids = []
        major_axis_lengths = []

        bb_areas = []
        convex_areas = []
        eccentricities = []
        equiv_diams = []
        extents = []
        minor_axis_lengths = []
        max_intensities = []
        mean_intensities = []
        min_intensities = []
        orientations = []
        perimeters = []
        solidities = []

        for pod in range(podocyte_count):

            areas.append(gen_props[pod].area)
            centroids.append(gen_props[pod].centroid)
            major_axis_lengths.append(gen_props[pod].major_axis_length)

            bb_areas.append(gen_props[pod].bbox_area)
            convex_areas.append(gen_props[pod].convex_area)
            eccentricities.append(gen_props[pod].eccentricity)
            equiv_diams.append(gen_props[pod].equivalent_diameter)
            extents.append(gen_props[pod].extent)
            minor_axis_lengths.append(gen_props[pod].minor_axis_length)
            max_intensities.append(gen_props[pod].max_intensity)
            mean_intensities.append(gen_props[pod].mean_intensity)
            min_intensities.append(gen_props[pod].min_intensity)
            orientations.append(gen_props[pod].orientation)
            perimeters.append(gen_props[pod].perimeter)
            solidities.append(gen_props[pod].solidity)

        pod_feat_vector.append(podocyte_count)
        area = np.mean(np.array(areas))*area_mpp2
        pod_feat_vector.append(area)
        bb_area = np.mean(np.array(bb_areas))*area_mpp2
        pod_feat_vector.append(bb_area)
        convex_area = np.mean(np.array(convex_areas))*area_mpp2
        pod_feat_vector.append(convex_area)
        eccentricity = np.mean(np.array(eccentricities))
        pod_feat_vector.append(eccentricity)
        equiv_diam = np.mean(np.array(equiv_diams))*dist_mpp
        pod_feat_vector.append(equiv_diam)
        extent = np.mean(np.array(extents))
        pod_feat_vector.append(extent)
        major_axis_length = np.mean(np.array(major_axis_lengths))*dist_mpp
        pod_feat_vector.append(major_axis_length)
        minor_axis_length = np.mean(np.array(minor_axis_lengths))*dist_mpp
        pod_feat_vector.append(minor_axis_length)
        max_intensity = np.mean(np.array(max_intensities))
        pod_feat_vector.append(max_intensity)
        mean_intensity = np.mean(np.array(mean_intensities))
        pod_feat_vector.append(mean_intensity)
        min_intensity = np.mean(np.array(min_intensities))
        pod_feat_vector.append(min_intensity)
        orientation = np.mean(np.array(orientations))
        pod_feat_vector.append(orientation)
        perimeter = np.mean(np.array(perimeters))*dist_mpp
        pod_feat_vector.append(perimeter)
        solidity = np.mean(np.array(solidities))
        pod_feat_vector.append(solidity)

        #write pod xml files
        #counter tool
        for pod in range(podocyte_count):
            centroid = np.array(centroids[pod]).reshape(-1,1)
#            print(centroid[0][0])
            xml_regionID = str(pod + 1)
            xml_Y = str(round(centroid[0][0]+x_start))
            xml_X = str(round(centroid[1][0]+y_start))

            xml_counter.write('\n\t\t\t<Region Id="' + xml_regionID + '" Type="5">\n\t\t\t\t<Vertices>\n\t\t\t\t\t<Vertex X="' + xml_X + '" Y="' + xml_Y + '" Z="0" />\n\t\t\t\t</Vertices>\n\t\t\t</Region>')

        #contours
        for pod in range(podocyte_count+1):
            if pod>0:
                pod_im = np.zeros(podocytes.shape)
                pod_im[podocyte_label==(pod)] = 1
                se = sk.morphology.disk(2)
                pod_im = sk.morphology.binary_dilation(pod_im,selem=se,out=None)
                pod_boundary = get_boundary(pod_im)

                L = []
                count = count+1
                xml_regionID = str(pod)
                index = pod-1
                length = major_axis_lengths[index]
                area = areas[index]
                length_um = length*dist_mpp
                area_um = area*area_mpp2
                xml_contour.write('\n\t\t\t<Region Id="' + str(count) + '" Type="0" Zoom="0" Selected="0" ImageLocation="" ImageFocus="-1" Length="' + str(length) + '" Area="'+ str(area) +'" LengthMicrons="'+ str(length_um) +'" AreaMicrons="'+ str(area_um) +'" Text="" NegativeROA="0" InputRegionId="0" Analyze="0" DisplayId="1">\n\t\t\t\t<Attributes/>\n\t\t\t\t<Vertices>\n')
                for point in pod_boundary:
                    xml_Y = str((point[1]+x_start))
                    xml_X = str((point[0]+y_start))
                    L.append(str('\t\t\t\t\t<Vertex X="' + xml_X + '" Y="' + xml_Y + '" Z="0"/>\n'))
                xml_contour.writelines(L)
                xml_contour.write('\t\t\t\t</Vertices>\n\t\t\t</Region>')

        #find bowmans and glom center distances
        glom_boundary = get_boundary(glom_mask)
        bowmans_dists = []
        glom_dists = []

        for pod in range(podocyte_count):
            centroid = np.array(centroids[pod]).reshape(1,-1)
            distances = sp.spatial.distance.cdist(centroid, glom_boundary, 'euclidean')
            bowmans_dists.append(np.min(distances))

            glom_label, glom_count = sp.ndimage.label(glom_mask)
            glom_centroid = sk.measure.regionprops(glom_label)
            glom_centroid = np.array(glom_centroid[0].centroid).reshape(1,-1)
            distances = sp.spatial.distance.cdist(centroid,glom_centroid, 'euclidean')
            glom_dists.append(np.min(distances))

        bowmans_dist = np.mean(np.array(bowmans_dists))*dist_mpp
        pod_feat_vector.append(bowmans_dist)
        glom_dist = np.mean(np.array(glom_dists))*dist_mpp
        pod_feat_vector.append(glom_dist)

        #find podocyte spatial density
        pod_dists = sp.spatial.distance.pdist(centroids)
        inter_pod_dist = sum(pod_dists)*dist_mpp
        if podocyte_count>1:
            pod_density = 1/inter_pod_dist
        else:
            pod_density = 0

        pod_feat_vector.append(inter_pod_dist)
        pod_feat_vector.append(pod_density)
        pod_feat_vector = np.array(pod_feat_vector)

    elif podocyte_count==0:
        pod_feat_vector = np.zeros([19,])

    return xml_counter, xml_contour, count, pod_feat_vector
