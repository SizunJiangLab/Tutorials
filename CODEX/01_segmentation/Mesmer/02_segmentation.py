import tifffile
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
from scipy.io import loadmat
import tensorflow as tf
from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay

'''
Date: 2024-08-27
Author: Huaying Qiu, Wenrui Wu
'''


# Setting ^ -------------------------------------------------------------------

# Constants
PIXEL_SIZE_UM = 0.5068

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Set up GPU configuration
os.environ["CUDA_VISIBLE_DEVICES"]="0"
try:
    tf_gpus = tf.config.list_physical_devices('GPU')
    for gpu in tf_gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
except:
    pass

# Setting $ -------------------------------------------------------------------


# Custom Function ^ -----------------------------------------------------------

def load_core_position(folder_input: str) -> pd.DataFrame:
    """
    Load core position data for a given TMA.
    Calculate the four vertices of each core.
    """
    path_dearrayer, = [str(path) for path in Path(folder_input).glob("*dearrayer*")]
    pos_df = pd.read_csv(path_dearrayer, sep='\t')
    
    diameter_mm, = re.findall(r"(?<=d=)[.0-9]+(?=[.]txt)", Path(path_dearrayer).name)
    diameter_mm = float(diameter_mm)
    radius_um = diameter_mm * 1000 / 2  # Convert diameter in mm to radius in um

    pos_df['x_beg_px'] = ((pos_df['Centroid X µm'] - radius_um) / PIXEL_SIZE_UM).apply(math.floor)
    pos_df['x_end_px'] = ((pos_df['Centroid X µm'] + radius_um) / PIXEL_SIZE_UM).apply(math.ceil)
    pos_df['y_beg_px'] = ((pos_df['Centroid Y µm'] - radius_um) / PIXEL_SIZE_UM).apply(math.floor)
    pos_df['y_end_px'] = ((pos_df['Centroid Y µm'] + radius_um) / PIXEL_SIZE_UM).apply(math.ceil)
    return pos_df

def rename_invalid_marker(marker: str):
    marker = re.sub(r'[/:]', '_', marker)
    return marker

def load_marker_list(folder_input: str) -> list:
    """
    Load interested marker lists.
    """
    im_name = Path(folder_input).name
    path_marker = f"{folder_input}/{im_name}_MarkerList.txt"

    with open(path_marker, 'r') as f:
        marker_list = f.read().splitlines()
    marker_list = [rename_invalid_marker(marker) for marker in marker_list if len(marker) > 0]
    return marker_list

def load_qptiff(folder_input: str) -> np.ndarray:
    """
    Load qptiff file for a given TMA.
    """
    im_name = Path(folder_input).name
    qptiff_path = f"{folder_input}/{im_name}.qptiff"
    qptiff_img = tifffile.imread(qptiff_path)
    return qptiff_img

def img_scale_marker(img: np.ndarray, marker: str, marker_list: List[str], scale: bool=True):
    img_marker = img[marker_list.index(marker)]
    if scale:
        img_marker = (img_marker - img_marker.min()) / (img_marker.max() - img_marker.min())
    return img_marker

def segment_mesmer_core(core_img: np.ndarray, 
                        boundary: List[str], internal: List[str], marker_list: List[str], 
                        scale: bool=True, 
                        maxima_threshold: float=0.075, interior_threshold: float=0.20):
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
        seg_stack, image_mpp = 0.5,
        postprocess_kwargs_whole_cell={
            'maxima_threshold': maxima_threshold,     # larger less cells
            'interior_threshold': interior_threshold  # larger larger cells
        }, 
        compartment = 'nuclear'
    )
    rgb_image = create_rgb_image(seg_stack, channel_colors = ["blue", "green"])
    overlay = make_outline_overlay(rgb_data = rgb_image, predictions = predictions)
    return predictions, overlay

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

def process_tma(folder_input: str, folder_output: str, 
                boundary: list, internal: list, 
                scale: bool=True, maxima_threshold: float=0.075, interior_threshold: float=0.20):
    """
    Process a single TMA.
    """
    img_name = Path(folder_input).name
    logging.info(f'Processing: {folder_input}')

    pos_df = load_core_position(folder_input)
    logging.info('Core positions loaded and vertices calculated')

    qptiff_img = load_qptiff(folder_input)
    logging.info('qptiff loaded for segmentation')

    marker_list = load_marker_list(folder_input)
    logging.info('Marker lists loaded')

    for i, row in tqdm(pos_df.iterrows()):
        core_name = row["Name"]
        core_img = qptiff_img[:, row["y_beg_px"]:row["y_end_px"], row["x_beg_px"]:row["x_end_px"]]

        folder_output_core = f"{folder_output}/{img_name}/{core_name}"
        os.makedirs(folder_output_core, exist_ok=True)
        
        # Segmentation
        segmentation_mask, overlay = segment_mesmer_core(core_img, boundary, internal, marker_list, scale, maxima_threshold, interior_threshold)
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

# Custom Function $ -----------------------------------------------------------


# Run Segmentation ^ ----------------------------------------------------------

def main():
    # select parameter ^ #####
    folder_input = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_raw/20240611_LMS-TMA_Scan1"
    folder_output = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_output/segmentation_run2/"
    # select parameter $ #####

    # select marker ^ #####
    boundary_1 = ['GLUT1', 'NaKATP', 'β-catenin', 'CD45', 'HLA1-300']
    boundary_2 = ['CD68', 'CD11c', 'Podoplanin', 'CD15', 'CD8', 'panCK', 'CD31', 'CD206', 'CD163']

    internal_1 = ['DAPI']
    internal_2 = ['α-SMA', 'Desmin', 'Vimentin']

    boundary = boundary_1 + boundary_2
    internal = internal_1 + internal_2
    # select marker $ #####

    # set parameter ^ #####
    scale = True
    maxima_threshold = 0.075
    interior_threshold = 0.20
    # set parameter $ #####

    # select folder ^ #####
    folder_input = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_raw/20240611_LMS-TMA_Scan1"
    folder_output = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_output/segmentation/"
    # select folder $ #####

    process_tma(
            folder_input, folder_output, 
            boundary, internal, 
            scale, maxima_threshold, interior_threshold
        )

if __name__ == "__main__":
    main()

# Run Segmentation $ ----------------------------------------------------------
