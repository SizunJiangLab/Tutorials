# Proseg Pipeline (Lite version)

https://github.com/dcjones/proseg

Proseg (**pro**babilistic **seg**mentation) is a cell segmentation method for in situ spatial transcriptomics. This pipeline is to help you set up the environment for CosMx segmentation using proseg. (This is a lite version of proseg pipeline without detailed explanation.)

**Notice**: Run commands within code chunks starting with `# Run in [somewhere] #####`. 


# Install Dependencies (only for the first time)

1. Rust programming language and Cargo (Rust package manager)
2. Julia programming language
3. Proseg
4. Executable shell script: `run_proseg.sh`

## 0. Set up working directory

```sh
# Run in shell #####

wd=/mnt/nfs/home/wenruiwu/Projects/proseg
mkdir -p $wd
```

## 1. Rust and Cargo

https://www.rust-lang.org/tools/install

```sh
# Run in shell #####

curl https://sh.rustup.rs -sSf | sh

# restart your current shell
```

## 2. Julia

https://github.com/JuliaLang/julia

```sh
# Run in shell #####

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

## 3. Proseg

https://github.com/dcjones/proseg

```sh
# Run in shell #####

# install proseg
cargo install proseg@1.1.7
```

## 4. `run_proseg.sh`

Download this shell script from our github repository `Tutorials/CosMx/01_segmentation/proseg/run_proseg.sh`.

```sh
# make the shell script executable
chmod a+x /path/to/run_proseg.sh
```


# Segmentation using proseg (executable shell script)

```sh
# Run in shell #####

# path of run_proseg.sh
run_proseg=/mnt/nfs/home/wenruiwu/Projects/proseg/output/run_proseg.sh

# parameter 1: directory for CosMx output files from AtoMx 
# parameter 2: directory for proseg output files
dir_cosmx=/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/
dir_output=/mnt/nfs/home/wenruiwu/Projects/proseg/output/TMA_542/

# run shell script
$run_proseg $dir_cosmx $dir_output
```

