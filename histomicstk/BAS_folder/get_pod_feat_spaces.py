"""
Created on Fri Sept 4 11:40:33 2020
@author: Briana Santo
"""
import numpy as np

def get_pod_feat_spaces(total_gloms):

    pod_feat_labels = ['num_pods','pod_areas','pod_bb_areas','pod_convex_areas','pod_eccentricities','pod_equiv_diams','pod_extents','pod_major_axis_lengths','pod_minor_axis_lengths','pod_max_intensities','pod_mean_intensities','pod_min_intensities','pod_orientations','pod_perimeters','pod_solidities','pod_bowmans_dists','pod_glom_dists','inter_pod_dists','pod_densities']
    pod_feat_qty = len(pod_feat_labels)
    pod_feat_array = np.empty([pod_feat_qty,total_gloms])

    return pod_feat_labels, pod_feat_qty, pod_feat_array
