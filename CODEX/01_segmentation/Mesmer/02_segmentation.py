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
import tensorflow as tf
import segmentation_mesmer


# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Set up GPU configuration
os.environ["CUDA_VISIBLE_DEVICES"]="3"
try:
    tf_gpus = tf.config.list_physical_devices('GPU')
    for gpu in tf_gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
except:
    pass

# Run Segmentation 
def main():
    # select parameter ^ #####
    folder_input = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_input/20240611_LMS-TMA_Scan1/"
    # select parameter $ #####

    segmentation_mesmer.process_tma(folder_input)

if __name__ == "__main__":
    main()