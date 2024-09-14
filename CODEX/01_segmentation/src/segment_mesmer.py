'''
Date: 2024-09-13
Author: Huaying Qiu, Wenrui Wu
'''

import tifffile
import json
import numpy as np
import pandas as pd 
from tqdm import tqdm 
import os 
from pathlib import Path
import logging
from typing import List, Dict
import skimage.io
import skimage.measure
import skimage.morphology
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay
from src.load_config import load_config

def load_tiff_markers(folder_marker: str) -> Dict[str, np.ndarray]:
    """
    Load all tiff of markers. 
    """
    marker_dict = {}
    for path_marker in Path(folder_marker).rglob('*.tiff'):
        marker_dict[str(path_marker.stem)] = tifffile.imread(path_marker)
    return marker_dict

def img_scale_marker(marker: str, marker_dict: Dict[str, np.ndarray], scale: bool=True) -> np.ndarray:
    """
    Scale image of a specific marker. 
    """
    img_marker = marker_dict[marker]
    if scale:
        img_marker = (img_marker - img_marker.min()) / (img_marker.max() - img_marker.min())
    return img_marker

def segment_mesmer(boundary: List[str], internal: List[str], marker_dict: Dict[str, np.ndarray], 
                   scale: bool=True, 
                   pixel_size_um: float=0.5068, maxima_threshold: float=0.075, interior_threshold: float=0.20):
    """
    Perform segmentation on a given image.
    """
    # Data for boundary markers
    boundary_list = [img_scale_marker(marker, marker_dict, scale=scale) for marker in boundary]
    boundary_sum = np.sum(boundary_list, axis=0)
    boundary_sum = 255 * (boundary_sum - boundary_sum.min()) / (boundary_sum.max() - boundary_sum.min())
    boundary_sum = boundary_sum.astype("uint8")
    # Data for internal markers
    internal_list = [img_scale_marker(marker, marker_dict, scale=scale) for marker in internal]
    internal_sum = np.sum(internal_list, axis=0)
    internal_sum = 255 * (internal_sum - internal_sum.min()) / (internal_sum.max() - internal_sum.min())
    internal_sum = internal_sum.astype("uint8")
    # Data for Mesmer
    seg_stack = np.stack((internal_sum, boundary_sum), axis = -1)
    seg_stack = np.expand_dims(seg_stack, 0)

    # Do whole cell segmentation
    # https://github.com/vanvalenlab/deepcell-tf/blob/master/notebooks/applications/Mesmer-Application.ipynb
    mesmer = Mesmer()
    segmentation_mask = mesmer.predict(
        seg_stack, image_mpp = pixel_size_um,
        postprocess_kwargs_whole_cell={
            'maxima_threshold': maxima_threshold,     # larger less cells
            'interior_threshold': interior_threshold  # larger larger cells
        }, 
        compartment = 'nuclear'
    )
    rgb_image = create_rgb_image(seg_stack, channel_colors = ["blue", "green"])
    overlay = make_outline_overlay(rgb_data = rgb_image, predictions = segmentation_mask)
    return segmentation_mask, rgb_image, overlay

def extract_single_cell_info(marker_dict: Dict[str, np.ndarray], 
                             segmentation_mask: np.ndarray):
    """
    Extract single cell information from a core.
    """
    marker_name = [marker for marker in marker_dict.keys() if marker != 'Empty']
    marker_array = [marker_dict[marker] for marker in marker_name]
    counts_no_noise = np.stack(marker_array, axis=2)

    stats = skimage.measure.regionprops(segmentation_mask)
    label_num = len(stats)

    channel_num = len(marker_array)
    data = np.zeros((label_num, channel_num))
    data_scale_size = np.zeros((label_num, channel_num))
    cell_sizes = np.zeros((label_num, 1))
    cell_props = np.zeros((label_num, 3))

    for i, region in enumerate(stats):
        cell_label = region.label
        label_counts = [counts_no_noise[coord[0], coord[1], :] for coord in region.coords]
        data[i] = np.sum(label_counts, axis=0)
        data_scale_size[i] = data[i] / region.area
        cell_sizes[i] = region.area
        cell_props[i] = [cell_label, region.centroid[0], region.centroid[1]]

    data_df = pd.DataFrame(data, columns=marker_name)
    data_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_df
    ], axis=1)

    data_scale_size_df = pd.DataFrame(data_scale_size, columns=marker_name)
    data_scale_size_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_scale_size_df
    ], axis=1)
    return data_full, data_scale_size_full

def segment_mesmer_object(folder_object: str, path_parameter: str, tag: str):
    """
    Process a single TMA.
    """
    name_object = Path(folder_object).name

    folder_output = f'{folder_object}/segmentation/{tag}/'
    os.makedirs(folder_output, exist_ok=True)    

    # Write parameter
    config = load_config(path_parameter)
    keys = ["boundary", "internal", "scale", "pixel_size_um",  "maxima_threshold", "interior_threshold"]
    config = {key: config.get(key) for key in keys}    
    with open(f'{folder_output}/parameter_segmentation.json', "w") as file:
        json.dump(config, file, indent=4)

    marker_dict = load_tiff_markers(f'{folder_object}/marker')
    logging.info('tiff loaded for segmentation')

    # Segmentation
    segmentation_mask, rgb_image, overlay = segment_mesmer(
        config["boundary"], config["internal"], marker_dict, 
        config["scale"], config["pixel_size_um"], config["maxima_threshold"], config["interior_threshold"]
    )
    segmentation_mask = segmentation_mask[0,...,0]
    overlay = overlay[0,...]
    tifffile.imwrite(f'{folder_output}/mesmer_mask.tiff', segmentation_mask)
    tifffile.imwrite(f'{folder_output}/mesmer_overlay.tiff', overlay)
    logging.info(f'Segmentation completed for {name_object}')
    
    # Extract single cell information
    data_full, data_scale_size_full = extract_single_cell_info(marker_dict, segmentation_mask)
    data_full.to_csv(f"{folder_output}/data.csv", index=False)
    data_scale_size_full.to_csv(f"{folder_output}/dataScaleSize.csv", index=False)
    logging.info(f'Single cell information extracted for {name_object}')