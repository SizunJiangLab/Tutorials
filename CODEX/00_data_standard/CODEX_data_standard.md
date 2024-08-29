# Data Standard for CODEX

**Date**: 2024-08-28

**Author**: Huaying Qiu, Wenrui Wu

---

## Summary

https://docs.google.com/presentation/d/1SGYQYQVoF5bTaQQvSDHbHqKvMe323egWCAPb7hHg9GY/edit#slide=id.g2f5dedbe25e_0_5

Collect everything of one TMA into one folder: 

- Folder: `path/to/folder/[name]`
    - Image: `[name].qptiff`
    - Marker list: `[name]_MarkerList.txt`
    - Core metadata: `[name]_metadata.csv`
    - TMA dearray information: `[name]_dearrayer.txt`
    - Parameter for segmentation: `[name]_parameter.json`

```
# Demo folder structure

/path/to/20240611_LMS-TMA_Scan1
├── 20240611_LMS-TMA_Scan1_dearrayer.txt
├── 20240611_LMS-TMA_Scan1_MarkerList.txt
├── 20240611_LMS-TMA_Scan1_metadata.csv
├── 20240611_LMS-TMA_Scan1_parameter.json
└── 20240611_LMS-TMA_Scan1.qptiff
```

## 1. Image

- Open the image using `QuPath` software. 
- Open `Image` tab and record the value of `Pixel width` or `Pixel height` (in μm, usually they are the same). 
- ⭐ Mark down the `Pixel width`/`Pixel height` into `[name]_parameter.json` file. 


## 2. Marker List

- Normally located at `/path/to/.temp/MarkerList.txt`. 
- Rename as `[name]_MarkerList.txt`. 


## 3. Core Metadata

- One row for one core. 
- One column for one category of information. 
- Include information of: 
    - `core_row`: row label, better to use alphabets (A, B, C, ...). 
    - `core_col`: column labels, better to use numbers (1, 2, 3, ...). 
    - `core_name`
    - `core_location`: should be `[core_row]-[core_col]` (A-1, B-2, C-3, ...).
    - Other information ...

### How to convert a wide table into a long table
- Input formula into one cell: `=ArrayFormula(FLATTEN(A2:A10&"|"&B1:H1&"|"&B2:H10))`
    - `A2:A10`: a range of cells for row labels. 
    - `B1:H1`: a range of cells for column labels. 
    - `B2:H10`: a range of cells for values (name or index of cores). 
    - `|`: separator. 
- Select the column with input formula. 
- `Data` > `Split text to columns` > `Separator` > `Custom` > `|`. 
- Add column names. 
- Input formula to calculate core location: `=CONCATENATE(J2,"-", K2)`
    - `J2`: cell location of `core_row`.
    - `K2`: cell location of `core_col`.
- Add other information. 
- Export the long table (`.csv` file) and rename as `[name]_metadata.csv`. 


## 4. TMA Dearray Information

- Open the image using `QuPath` software. 
- `TMA` > `TMA dearrayer`: 
    - `TMA core diameter`: if you don’t know, you can adjust the value base on the result. 
    - `Row labels` and `Column labels`: should be the same as those in metadata (better to use alphabets for row labels and numbers for column labels). 
    - `Bounds scale factor`: set to `100`. 
    - Click `Run`. 
- Adjust location of the circles:
    - Double click to select a not-well-aligned circle (the selected circle will turn yellow). 
    - Drag it and move to an appropriate location. 
    - Check every core to make sure the circle contains your entire core. 
- Add missing rows or columns: 
    - Double click to select a circle above or below the missing row or column. 
    - Right click the selected circle > `TMA` > `Add` > `Add TMA row/column before/after`. (After you add/delete a row/column, you are asked to relabel the TMA grid. Just use the labels at the beginning. )
    - Double click to select a newly added circle (which is hidden when added), right click it > `TMA` > `Set core valid`. 
    - Repeat for all missing circles and move them to appropriate locations. 
- `Measure` > `Show TMA measurements` > `Save` as `[name]_dearrayer.txt`. 
- ⭐ Mark down the value of `TMA core diameter` into `[name]_parameter.json` file. 
- (Optional) In case you need to modify anything in the future and you don't want to align everything again, you can save your TMA grid as a `.qpdata` file. 
    - `File` > `Save` or `Save as`. 
    - Next time you want to reopen it, just open the `.qpdata` file using `QuPath`. 

## 5. Parameter for segmentation

- `folder_input`: path of the folder following this data standard, which contains all files needed for segmentation. 
- `folder_output`: path of the folder to store output files of segmentation. Folder named as `[name]` for each TMA and sub-folders named as `[core_name]` for each core will be created. 
```
/path/to/folder_output
├── 20240611_LMS-TMA_Scan1
│   ├── A-1
│   ├── B-1
│   └── ...
├── 20240611_TMA_SM1-4_Scan1
│   ├── A-1
│   ├── B-1
│   └── ...
└── ...
```
- `diameter_mm`: diameter of dearrayer circles in millimeter. 
- `boundary`: names of markers as "boundary markers" used for segmentation, which are usually membrane markers. 
- `internal`: names of markers as "internal markers" used for segmentation, which are usually nucleus and cytoplasmic markers.
- `scale`: `true` to scale the signals of `boundary` and `internal` markers before summing up. 
- `pixel_size_um`: pixel size of `.qptiff` image in micrometer.  
- `maxima_threshold`: parameter for segmentation, larger value for less cells. 
- `interior_threshold`: parameter for segmentation, larger value for larger cells.   


```json
# Demo
{
    "folder_input": "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_input/20240611_LMS-TMA_Scan1",
    "folder_output": "/mnt/nfs/home/wenruiwu/projects/Sarcoma/data_test/data_output", 

    "diameter_mm": 1.3,
    
    "boundary": ["GLUT1", "NaKATP", "β-catenin", "CD45", "HLA1-300", "CD68", "CD11c", "Podoplanin", "CD15", "CD8", "panCK", "CD31", "CD206", "CD163"], 
    "internal": ["DAPI", "α-SMA", "Desmin", "Vimentin"], 
    "scale": true, 
    "pixel_size_um": 0.5068,
    "maxima_threshold": 0.075, 
    "interior_threshold": 0.20
}
```


