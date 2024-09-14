'''
Date: 2024-09-13
Author: Huaying Qiu, Wenrui Wu
'''


import tifffile
import json
import numpy as np
import math 
import pandas as pd 
from tqdm import tqdm 
import os 
import re
import logging
from typing import List, Dict
from src.load_config import load_config


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

def crop_tma(path_parameter):
    config = load_config(path_parameter)
    keys = ["name", "path_qptiff", "path_dearrayer", "path_marker",  "pixel_size_um", "diameter_mm", "folder_output"]
    config = {key: config.get(key) for key in keys}
    with open(f'{config["folder_output"]}/{config["name"]}/parameter_crop.json', "w") as file:
        json.dump(config, file, indent=4)

    pos_df = load_core_position(config["path_dearrayer"], config["pixel_size_um"], config["diameter_mm"])
    logging.info('Core positions loaded and vertices calculated')
    qptiff_img = load_qptiff(config["path_qptiff"])
    logging.info('qptiff loaded for segmentation')
    marker_list = load_marker_list(config["path_marker"])
    logging.info('Marker lists loaded')

    for _, row in tqdm(pos_df.iterrows()):
        core_name = row["Name"]
        core_img = qptiff_img[:, row["y_beg_px"]:row["y_end_px"], row["x_beg_px"]:row["x_end_px"]]
        folder_output_core = f'{config["folder_output"]}/{config["name"]}/{core_name}/marker'
        os.makedirs(folder_output_core, exist_ok=True)
        # Save markers
        for marker in marker_list:
            tifffile.imwrite(f"{folder_output_core}/{marker}.tiff", core_img[marker_list.index(marker)])
        logging.info(f'Markers saved for core {core_name}')