# Proseg Pipline

https://github.com/dcjones/proseg

- Proseg (**pro**babilistic **seg**mentation) is a cell segmentation method for in situ spatial transcriptomics.
- This pipeline is to help you set up the environment for CosMx segmentation using proseg. 


# Install Dependencies (only for the first time)

1. Rust programming language and Cargo (Rust package manager)
2. Julia programming language
3. proseg

## 0. Set up working directory

```sh
wd=/mnt/nfs/home/wenruiwu/Projects/proseg
mkdir $wd
```

## 1. Rust and Cargo

https://www.rust-lang.org/tools/install

```sh
curl https://sh.rustup.rs -sSf | sh

# restart your current shell
```

## 2. Julia

https://github.com/JuliaLang/julia

```sh
# clone Julia repository
cd $wd
git clone https://github.com/JuliaLang/julia.git
cd $wd/julia

# use the most recent stable version (check the julia repository)
git checkout v1.10.4

# build the julia executable
make

# add Julia to your $PATH (only run it at the first time)
echo export PATH="$PATH:$wd/julia" >> ~/.bashrc

# restart your current shell

# install packages in Julia
julia -e 'import Pkg; Pkg.add(["Glob", "CSV", "DataFrames", "CodecZlib", "ArgParse"])'
```

## 3. proseg

https://github.com/dcjones/proseg

```sh
# install proseg
cargo install proseg

# clone proseg repository
cd $wd
git clone https://github.com/dcjones/proseg.git
```


# Segmentation using proseg (original)

## 1. Fix bugs regarding input files in `stitch-cosmx.jl`

(2024-07-28) Codes regarding the file structure of CosMx output in the proseg repository is incorrect (due to an update of AtoMx changing the output file structure). There is already a issue in the proseg repository (https://github.com/dcjones/proseg/issues/3) planning to fix it. Before the author release an updated version fixing this bug, we need to fix it by ourselves. 

```sh
# diretory of proseg repository
dir_proseg=/mnt/nfs/home/wenruiwu/Projects/proseg/proseg

sed 's/glob("S0\/\*\//glob("/g' \
    $dir_proseg/extra/stitch-cosmx.jl > \
    $dir_proseg/extra/stitch-cosmx_modified.jl
```

## 2. Set path of input and output

```sh
# path of stitch-cosmx_modified.jl
stitch_cosmx=/mnt/nfs/home/wenruiwu/Projects/proseg/proseg/extra/stitch-cosmx_modified.jl

# directory for raw CosMx files 
# should be direct parent directory of RunSummary and AnalysisResults folders
dir_rawfile=/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/RawFiles/RCC_BPC_23_tma542/20240709_223954_S1

# directory for proseg output
dir_output=/mnt/nfs/home/wenruiwu/Projects/proseg/output/TMA_542
```

## 3. Run porseg

`--voxel-layers 5`: To use 5 layers of voxels on the z-axis. Essentially how 3D the segmentation should be. (The thickness of our slides for CosMx is 5 μm. If we separate the z-axis into 5 layer, the length of voxel on z-axis is 1 μm. By default, the length of voxel on x-axis and y-axis is 1 μm. So the volume of each voxel is 1 μm<sup>3</sup>.) 

```sh
mkdir $dir_output
cd $dir_output
julia $stitch_cosmx $dir_rawfile transcripts.csv.gz
proseg --cosmx-micron --voxel-layers 5 transcripts.csv.gz
```

# Segmentation using proseg (executable shell script)

## 1. Introduction

Because the output file structure of AtoMx is messy, it is easy to input the long and messy diretory incorrectly. Morevoer, we need to fix the bugs regarding input files in `stitch-cosmx.jl`. So I create an executable shell script to make the proseg convenient to run. 

The file structure of the AtoMx outfile is like:

```
├── /mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132 
│   ├── RawFiles 
│   │   └── RCC_BPC_23_tma542 
│   │       ├── 20240709_223954_S1 
│   │       │   ├── RunSummary 
│   │       │   ├── AnalysisResults 
│   │       │   └── ...
│   │       └── Logs 
│   │           └── ...
│   ├── flatFiles
│   │   └── RCC_BPC_23_tma542 
│   │       └── 20240709_223954_S1 
│   │           ├── RCC_BPC_23_tma542_exprMat_file.csv.gz*
│   │           ├── RCC_BPC_23_tma542_fov_positions_file.csv.gz*
│   │           ├── RCC_BPC_23_tma542_metadata_file.csv.gz*
│   │           ├── RCC_BPC_23_tma542-polygons.csv.gz*
│   │           └── RCC_BPC_23_tma542_tx_file.csv.gz*
│   └── ... 
```

- In the *original* method, the input directory should be **direct parent directory** of RunSummary and AnalysisResults folders: `/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/RawFiles/RCC_BPC_23_tma542/20240709_223954_S1`
- While in our *executable shell script* method, we just input the directory containing the RunSummary and AnalysisResults folders: `/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132`. The script will automatically search the direct parent directory of them. 

## 2. Run `run_proseg.sh`

```sh
# path of run_proseg.sh
run_proseg=/mnt/nfs/home/wenruiwu/Projects/proseg/output/run_proseg.sh

# make the shell script executable
chmod a+x $run_proseg

# run shell script
## parameter 1: path of stitch-cosmx.jl
## parameter 2: directory for CosMx output files from AtoMx 
## parameter 3: directory for proseg output files
$run_proseg \
    /mnt/nfs/home/wenruiwu/Projects/proseg/proseg/extra/stitch-cosmx.jl \
    /mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/ \
    /mnt/nfs/home/wenruiwu/Projects/proseg/output/TMA_542/
```

