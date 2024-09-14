'''
Date: 2024-09-12
Author: Huaying Qiu, Wenrui Wu
'''

import sys
import logging
from src.crop_tma import crop_tma

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    if len(sys.argv) != 2:
        print("Usage: python 02_crop_tma.py <path_parameter>")
        sys.exit(1)
    path_parameter = sys.argv[1]
    # path_parameter = "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_input/20240611_LMS-TMA_Scan1/20240611_LMS-TMA_Scan1_parameter.json"
    crop_tma(path_parameter)

if __name__ == "__main__":
    main()
