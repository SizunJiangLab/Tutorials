'''
Date: 2024-08-28
Author: Huaying Qiu, Wenrui Wu
'''


import tifffile
import json
import numpy as np
import math 
import pandas as pd 
from tqdm import tqdm 
import os 
from pathlib import Path
import re
import logging
from typing import List, Dict
import skimage.io
import skimage.measure
import skimage.morphology
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay


def get_path_file(folder_input: str) -> dict:
    """
    Get pathes of all files needed within the input folder. 
    """
    name = Path(folder_input).name
    path_dict = {
        "name": name, 
        "path_qptiff": f"{folder_input}/{name}.qptiff", 
        "path_marker": f"{folder_input}/{name}_MarkerList.txt", 
        "path_metadata": f"{folder_input}/{name}_metadata.txt", 
        "path_dearrayer": f"{folder_input}/{name}_dearrayer.txt"
    }
    return path_dict
    
def load_config(folder_input: str) -> dict:
    """
    Load all the paramters needed. 
    """
    img_name = Path(folder_input).name
    path_config = f"{folder_input}/{img_name}_parameter.json"
    with open(path_config, 'r') as f:
        config_dict = json.load(f)

    get_all_path = config_dict["folder_input"]
    path_dict = get_path_file(get_all_path)
    return config_dict | path_dict

def load_core_position(path_dearrayer: str, pixel_size_um: float, diameter_mm: float) -> pd.DataFrame:
    """
    Load core position data for a given TMA.
    Calculate the four vertices of each core.
    """
    pos_df = pd.read_csv(path_dearrayer, sep='\t')
    radius_um = diameter_mm * 1000 / 2  # Convert diameter in mm to radius in um

    pos_df['x_beg_px'] = ((pos_df['Centroid X µm'] - radius_um) / pixel_size_um).apply(math.floor)
    pos_df['x_end_px'] = ((pos_df['Centroid X µm'] + radius_um) / pixel_size_um).apply(math.ceil)
    pos_df['y_beg_px'] = ((pos_df['Centroid Y µm'] - radius_um) / pixel_size_um).apply(math.floor)
    pos_df['y_end_px'] = ((pos_df['Centroid Y µm'] + radius_um) / pixel_size_um).apply(math.ceil)
    return pos_df

def rename_invalid_marker(marker: str) -> str:
    """
    rename the invalid marker name (containing "/" and ":"). 
    """
    marker = re.sub(r'[/:]', '_', marker)
    return marker

def load_marker_list(path_marker: str) -> list:
    """
    Load interested marker lists.
    """
    with open(path_marker, 'r') as f:
        marker_list = f.read().splitlines()
    marker_list = [rename_invalid_marker(marker) for marker in marker_list if len(marker) > 0]
    return marker_list

def load_qptiff(path_qptiff: str) -> np.ndarray:
    """
    Load qptiff file for a given TMA.
    """
    qptiff_img = tifffile.imread(path_qptiff)
    return qptiff_img

def img_scale_marker(img: np.ndarray, marker: str, marker_list: List[str], scale: bool=True) -> np.ndarray:
    img_marker = img[marker_list.index(marker)]
    if scale:
        img_marker = (img_marker - img_marker.min()) / (img_marker.max() - img_marker.min())
    return img_marker

def segment_mesmer_core(core_img: np.ndarray, 
                        boundary: List[str], internal: List[str], marker_list: List[str], 
                        scale: bool=True, 
                        pixel_size_um: float=0.5068, maxima_threshold: float=0.075, interior_threshold: float=0.20):
    """
    Perform segmentation on a single core image.
    """
    # Prepare data
    boundary = [rename_invalid_marker(marker) for marker in boundary]
    boundary_list = [img_scale_marker(core_img, marker, marker_list, scale=scale) for marker in boundary]
    boundary_sum = np.sum(boundary_list, axis=0)
    boundary_sum = 255 * (boundary_sum - boundary_sum.min()) / (boundary_sum.max() - boundary_sum.min())
    boundary_sum = boundary_sum.astype("uint8")

    internal = [rename_invalid_marker(marker) for marker in internal]
    internal_list = [img_scale_marker(core_img, marker, marker_list, scale=scale) for marker in internal]
    internal_sum = np.sum(internal_list, axis=0)
    internal_sum = 255 * (internal_sum - internal_sum.min()) / (internal_sum.max() - internal_sum.min())
    internal_sum = internal_sum.astype("uint8")

    seg_stack = np.stack((internal_sum, boundary_sum), axis = -1)
    seg_stack = np.expand_dims(seg_stack, 0)

    # Do whole cell segmentation
    # https://github.com/vanvalenlab/deepcell-tf/blob/master/notebooks/applications/Mesmer-Application.ipynb
    mesmer = Mesmer()
    predictions = mesmer.predict(
        seg_stack, image_mpp = pixel_size_um,
        postprocess_kwargs_whole_cell={
            'maxima_threshold': maxima_threshold,     # larger less cells
            'interior_threshold': interior_threshold  # larger larger cells
        }, 
        compartment = 'nuclear'
    )
    rgb_image = create_rgb_image(seg_stack, channel_colors = ["blue", "green"])
    overlay = make_outline_overlay(rgb_data = rgb_image, predictions = predictions)
    return predictions, rgb_image, overlay

def extract_single_cell_info(core_dict: Dict[str, np.ndarray], 
                             segmentation_mask: np.ndarray, 
                             marker_list: List[str]):
    """
    Extract single cell information from a core.
    """
    array_list = [core_dict[channel] for channel in marker_list if channel != 'Empty']
    counts_no_noise = np.stack(array_list, axis=2)

    stats = skimage.measure.regionprops(segmentation_mask)
    label_num = len(stats)

    channel_num = len(array_list)
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

    col_names = [marker for marker in marker_list if marker != 'Empty']

    data_df = pd.DataFrame(data, columns=col_names)
    data_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_df
    ], axis=1)

    data_scale_size_df = pd.DataFrame(data_scale_size, columns=col_names)
    data_scale_size_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_scale_size_df
    ], axis=1)
    return data_full, data_scale_size_full

def process_tma(folder_input: str):
    """
    Process a single TMA.
    """
    logging.info(f'Processing: {folder_input}')
    config = load_config(folder_input)

    pos_df = load_core_position(config["path_dearrayer"], config["pixel_size_um"], config["diameter_mm"])
    logging.info('Core positions loaded and vertices calculated')

    qptiff_img = load_qptiff(config["path_qptiff"])
    logging.info('qptiff loaded for segmentation')

    marker_list = load_marker_list(config["path_marker"])
    logging.info('Marker lists loaded')

    for i, row in tqdm(pos_df.iterrows()):
        core_name = row["Name"]
        core_img = qptiff_img[:, row["y_beg_px"]:row["y_end_px"], row["x_beg_px"]:row["x_end_px"]]

        folder_output_core = f"{config['folder_output']}/{config['name']}/{core_name}"
        os.makedirs(folder_output_core, exist_ok=True)
        
        # Segmentation
        segmentation_mask, rgb_image, overlay = segment_mesmer_core(
            core_img, 
            config["boundary"], config["internal"], marker_list, 
            config["scale"], config["pixel_size_um"], config["maxima_threshold"], config["interior_threshold"]
        )
        segmentation_mask = segmentation_mask[0,...,0]
        overlay = overlay[0,...]
        tifffile.imwrite(f"{folder_output_core}/mesmer_mask.tiff", segmentation_mask)
        tifffile.imwrite(f"{folder_output_core}/mesmer_overlay.tiff", overlay)
        logging.info(f'Segmentation completed for core {core_name}')

        # Save markers
        core_dict = {}
        for marker in marker_list:
            tifffile.imwrite(f"{folder_output_core}/{marker}.tiff", core_img[marker_list.index(marker)])
            core_dict[marker] = core_img[marker_list.index(marker)]   
        logging.info(f'Markers saved for core {core_name}')

        # Extract single cell information
        data_full, data_scale_size_full = extract_single_cell_info(core_dict, segmentation_mask, marker_list)
        data_full.to_csv(f"{folder_output_core}/data.csv", index=False)
        data_scale_size_full.to_csv(f"{folder_output_core}/dataScaleSize.csv", index=False)
        logging.info(f'Single cell information extracted for core {core_name}')

# [todo] read and cut mask/marker (if shape > 3000 * 3000)
# def cut_and_save_cores(tma: str, core_dict: Dict):
#     """
#     Cut cores into smaller sections and save as TIFF files.
#     """
#     for core_num, markers in tqdm(core_dict.items(), desc="Cutting and saving cores"):
#         seg_path = f'{OUTPUT_DIR}/seg_results_8bit/{tma}/{core_num}/MESMER_mask.tiff'
#         try:
#             seg_mask = tifffile.imread(seg_path)
#         except FileNotFoundError:
#             logging.warning(f"Segmentation mask not found for core {core_num}. Skipping.")
#             continue
#
#         img_height, img_width = seg_mask.shape
#         height_mid, width_mid = img_height // 2, img_width // 2
#
#         for fov in range(4):
#             output_path = f'{OUTPUT_DIR}/mantis_img_8bit/{tma}/{core_num}_fov{fov}'
#             os.makedirs(output_path, exist_ok=True)
#
#             x1 = (fov % 2) * width_mid
#             x2 = ((fov + 1) % 2) * width_mid + ((fov + 2) % 2) * img_width
#             y1 = (fov // 2) * height_mid
#             y2 = ((fov // 2 + 1) % 2) * height_mid + ((fov // 2 + 2) % 2) * img_height
#
#             tifffile.imwrite(os.path.join(output_path, 'seg_mask.tiff'), seg_mask[y1:y2, x1:x2])
#
#             for marker, img in markers.items():
#                 marker_fov_img = img[y1:y2, x1:x2]
#                 tifffile.imwrite(os.path.join(output_path, f'{marker}.tiff'), marker_fov_img)