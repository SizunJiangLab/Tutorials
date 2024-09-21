'''
Date: 2024-09-21
Author: Huaying Qiu, Wenrui Wu
'''

import json
from pathlib import Path

def load_config(path_parameter: str) -> dict:
    """
    Load all the parameters needed. 
    """
    with open(path_parameter, "r") as f:
        config_dict = json.load(f)
    dir_input = config_dict["dir_input"]
    dir_output = config_dict["dir_output"]
    name = config_dict["name"]
    path_dict = {
        "name": name, 
        "folder_input": Path(dir_input, name), 
        "folder_output": Path(dir_output, name), 
        "path_qptiff": Path(dir_input, name, f"{name}.qptiff"), 
        "path_marker": Path(dir_input, name, f"{name}_MarkerList.txt"), 
        "path_metadata": Path(dir_input, name, f"{name}_metadata.csv"), 
        "path_dearrayer": Path(dir_input, name, f"{name}_dearrayer.txt"), 
    }
    path_dict = {key: str(value) for key, value in path_dict.items()}
    return config_dict | path_dict