# directory for raw CosMx files 
dir_cosmx=/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/
# directory for proseg output files
dir_output=/mnt/nfs/home/wenruiwu/Projects/proseg/output/TMA_542/

dir_cosmx=$1
dir_output=$2

mkdir -p $dir_output
cd $dir_output

git clone https://github.com/dcjones/proseg.git
stitch_cosmx=$dir_output/proseg/extra/stitch-cosmx.jl
# edit the file paths in stitch-cosmx.jl file
sed 's/glob("S0\/\*\//glob("/g' \
    $stitch_cosmx > \
    $dir_output/stitch-cosmx_modified.jl
stitch_cosmx=$dir_output/stitch-cosmx_modified.jl
rm -rf proseg

# find the parent directory of RunSummary and AnalysisResults folders
dir_cosmx=$(dirname $(find $dir_cosmx -name RunSummary -type d -print -quit))

julia $stitch_cosmx $dir_cosmx transcripts.csv.gz
proseg --cosmx-micron --voxel-layers 5 transcripts.csv.gz