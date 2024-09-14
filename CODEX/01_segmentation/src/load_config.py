'''
Date: 2024-09-13
Author: Huaying Qiu, Wenrui Wu
'''


import json
from pathlib import Path


def get_path_file(folder_input: str) -> dict:
    """
    Get paths of all files needed within the input folder. 
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
    
def load_config(path_parameter: str) -> dict:
    """
    Load all the parameters needed. 
    """
    with open(path_parameter, "r") as f:
        config_dict = json.load(f)
    get_all_path = config_dict["folder_input"]
    path_dict = get_path_file(get_all_path)
    return config_dict | path_dict