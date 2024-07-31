# Proseg Pipeline (Notes)

https://github.com/dcjones/proseg

Some notes for detailed information of the proseg pipeline. 


# Bug in `stitch-cosmx.jl` (2024-07-28) 

## 1. What is the bug? 

Codes regarding the file structure of CosMx output in the proseg repository is incorrect (due to an update of AtoMx changing the output file structure). There is already a issue in the proseg repository (https://github.com/dcjones/proseg/issues/3) planning to fix it. Before the author release an updated version fixing this bug, we need to fix it by ourselves. 

The latest file structure of the AtoMx output file is like:

```sh
# file with (*) is what proseg needs
├── /mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132 
│   ├── RawFiles 
│   │   └── RCC_BPC_23_tma542 
│   │       ├── 20240709_223954_S1 
│   │       │   └── RunSummary 
│   │       │   │   ├── (*) *_ExptConfig.txt
│   │       │   │   ├── (*) latest.fovs.csv
│   │       │   │   └── ...
│   │       │   ├── AnalysisResults 
│   │       │   │   └── xvjom9x9tx
│   │       │   │       ├── (*) FOV00001
│   │       │   │       ├── (*) FOV00002
│   │       │   │       └── ...
│   │       │   └── ...
│   │       └── Logs 
│   │           └── ...
│   ├── flatFiles
│   │   └── RCC_BPC_23_tma542 
│   │       └── 20240709_223954_S1 
│   │           ├── RCC_BPC_23_tma542_exprMat_file.csv.gz
│   │           ├── RCC_BPC_23_tma542_fov_positions_file.csv.gz
│   │           ├── RCC_BPC_23_tma542_metadata_file.csv.gz
│   │           ├── RCC_BPC_23_tma542-polygons.csv.gz
│   │           └── RCC_BPC_23_tma542_tx_file.csv.gz
│   └── ... 
```

In the `stitch-cosmx.jl`, the proseg uses regular expression to search the files it needs:

```julia
config_filename = glob("S0/*/RunSummary/*_ExptConfig.txt", path)[1]
fov_filename = glob("S0/*/RunSummary/latest.fovs.csv", path)[1]
fov_paths = glob("S0/*/AnalysisResults/*/FOV*", path)
```

The regular expressions indicate that the file structure is as followed, which has changed after the update of AtoMx:  

```sh
# file with (*) is what proseg needs
├── /mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132 
│   ├── S0 
│   │   └── * 
│   │       ├── RunSummary
│   │       │   ├── (*) *_ExptConfig.txt
│   │       │   ├── (*) latest.fovs.csv
│   │       │   └── ...
│   │       ├── AnalysisResults 
│   │       │   │   └── xvjom9x9tx
│   │       │   │       ├── (*) FOV00001
│   │       │   │       ├── (*) FOV00002
│   │       │   │       └── ...
│   │       │   └── ...
│   │       └── ... 
│   └── ... 
```

As the file structure has been changed, the regular expressions in `stitch-cosmx.jl` cannot get the correct path for the files we need. 


## 2. How to fix this bug?

To match the correct files that the proseg needs, we create a modified script named `stitch-cosmx_modified.jl`. We replace  `glob("S0/*/` with `glob("` in `stitch-cosmx.jl` (It is the same as just removing the `S0/*/`, but we use `glob("` to locate the specific regular expression we want to modify.)

```sh
# Run in shell #####

# directory of proseg repository
dir_proseg=/mnt/nfs/home/wenruiwu/Projects/proseg/proseg

sed 's/glob("S0\/\*\//glob("/g' \
    $dir_proseg/extra/stitch-cosmx.jl > \
    $dir_proseg/extra/stitch-cosmx_modified.jl
```

With the modified regular expressions, we can match the correct files that the proseg needs: 

```julia
# path: direct parent directory of RunSummary and AnalysisResults folders
config_filename = glob("RunSummary/*_ExptConfig.txt", path)[1]
fov_filename = glob("RunSummary/latest.fovs.csv", path)[1]
fov_paths = glob("AnalysisResults/*/FOV*", path)
```

# Segmentation using proseg (executable shell script)

## Some inconveniences of the original proseg pipeline
 
 - We need to fix the bug regarding input files in `stitch-cosmx.jl`.
 - The output file structure of AtoMx is messy, it is easy to input the long and messy directory incorrectly. 

## Conveniences of executable shell script

- Automatically fix the bug in `stitch-cosmx.jl`
- Easy input of parameters

    - In the *original* method, the input directory should be **direct parent directory** of RunSummary and AnalysisResults folders: `/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/RawFiles/RCC_BPC_23_tma542/20240709_223954_S1`
    - While in our *executable shell script* method, we just input the directory containing the RunSummary and AnalysisResults folders: `/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132`. The script will automatically search the direct parent directory of them. 