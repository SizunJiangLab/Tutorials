#!/bin/bash

# edit your parameter ^ #####
path_parameter=/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_input/20240611_LMS-TMA_Scan1/20240611_LMS-TMA_Scan1_parameter.json
path_objects=/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_output/20240611_LMS-TMA_Scan1
tag=test_20240914
# edit your parameter & #####

# Activate cellSeg environment 
source /opt/miniconda3/bin/activate cellSeg

# Crop TMA
python \
    /mnt/nfs/home/wenruiwu/projects/Tutorials/CODEX/01_segmentation/02_crop_tma.py \
    ${path_parameter}

# Segmentation
python /mnt/nfs/home/wenruiwu/projects/Tutorials/CODEX/01_segmentation/03_segment_mesmer_tma.py \
    ${path_objects} \
    ${path_parameter} \
    ${tag}
