'''
Date: 2024-08-28
Author: Huaying Qiu, Wenrui Wu
'''

import sys
import os 
import logging
import tensorflow as tf
from src.segment_mesmer import segment_mesmer_object

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Set up GPU configuration
os.environ["CUDA_VISIBLE_DEVICES"]="3"
try:
    tf_gpus = tf.config.list_physical_devices('GPU')
    for gpu in tf_gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    logging.info("Succeed to set up GPU")
except:
    logging.info("Fail to set up GPU")
    pass

# Run Segmentation 
def main():
    if len(sys.argv) != 4:
        print("Usage: python 03_segment_mesmer.py <folder_object>, <path_parameter>, <tag>")
        sys.exit(1)
    folder_object = sys.argv[1]
    path_parameter = sys.argv[2]
    tag = sys.argv[3]
    segment_mesmer_object(folder_object, path_parameter, tag)

if __name__ == "__main__":
    main()
