# path of stitch-cosmx.jl
stitch_cosmx=/mnt/nfs/home/wenruiwu/Projects/proseg/proseg/extra/stitch-cosmx.jl
# directory for raw CosMx files 
dir_rawfile=/mnt/nfs/storage/CosMX/RCC_TMA542_section5_v132/
# directory for proseg output files
dir_output=/mnt/nfs/home/wenruiwu/Projects/proseg/output/TMA_542/

stitch_cosmx=$1
dir_rawfile=$2
dir_output=$3

# edit the file paths in stitch-cosmx.jl file
sed 's/glob("S0\/\*\//glob("/g' \
    $stitch_cosmx > \
    $(dirname $stitch_cosmx)/stitch-cosmx_modified.jl
stitch_cosmx=$(dirname $stitch_cosmx)/stitch-cosmx_modified.jl

# find the parent directory of RunSummary and AnalysisResults folders
path_rawfile=$(dirname $(find $dir_rawfile -name RunSummary -type d -print -quit))

mkdir -p $dir_output
cd $dir_output
julia $stitch_cosmx $path_rawfile transcripts.csv.gz
proseg --cosmx-micron --voxel-layers 5 transcripts.csv.gz