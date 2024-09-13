'''
Date: 2024-09-12
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

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

path_parameter = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_input/20240611_LMS-TMA_Scan1/20240611_LMS-TMA_Scan1_parameter.json"
segmentation_mesmer.crop_tma(path_parameter) 