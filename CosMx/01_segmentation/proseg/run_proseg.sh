#!/bin/bash

# Script to process CosMx data using proseg
#
# Usage:
#   ./run_proseg.sh [dir_cosmx] [dir_output] [--voxel-layers n_layer] [--nthreads n_threads]
#
# Arguments:
#   dir_cosmx       Directory containing CosMx output exported from AtoMx SIP, including
#                   "*_tx_file.csv.gz" file.
#   dir_output      Directory to store proseg output files.
#   --voxel-layers  (Optional) Set the number of voxel layers manually. If not provided,
#                   it will be determined automatically. In automatic mode, the script
#                   reads the $tx_file to find the 'z' column, calculates the unique values
#                   in that column, and sets the number of voxel layers to match the number
#                   of unique 'z' values.
#   --nthreads       (Optional) Set the number of threads to use. Default is 28.
#
# Description:
#   This script processes CosMx transcript data using the proseg tool. It takes in a
#   directory containing CosMx output files and an output directory where results will be
#   stored. The script will search for the required transcript file, and if found, it will
#   run proseg with the necessary parameters to generate various output files, including
#   expected counts, cell metadata, transcript metadata, and cell polygons. If the
#   --voxel-layers argument is not provided, the script will determine the appropriate
#   number of layers by analyzing the 'z' column in the transcript file, ensuring the
#   correct number of layers based on the data's z-depth.

# Verify if the dir_cosmx argument is provided
if [ -z "$1" ]; then
    echo "Error: Please provide the dir_cosmx as the first argument."
    exit 1
else
    dir_cosmx="$1"
fi

# Verify if the dir_output argument is provided
if [ -z "$2" ]; then
    echo "Error: Please provide the dir_output as the second argument."
    exit 1
else
    dir_output="$2"
fi

# Set default number of threads
n_threads=28

# Parse optional arguments
n_layer="auto"
shift 2
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --voxel-layers)
            n_layer="$2"
            shift 2
            ;;
        --nthreads)
            n_threads="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create the output directory if it doesn't exist
mkdir -p "$dir_output"
cd "$dir_output"

# Find the specific file based on the pattern
tx_file=$(find "$dir_cosmx" -type f -name "*_tx_file.csv.gz")

# Check if the file was found
if [ -z "$tx_file" ]; then
    echo "Error: No matching file found in $dir_cosmx."
    exit 1
else
    echo "File found: $tx_file"
fi

# Determine voxel layers if set to auto
if [ "$n_layer" == "auto" ]; then
    echo "Determining number of voxel layers automatically..."
    # Determine the column index for 'z'
    z_column=$(zcat "$tx_file" | head -n 1 | tr ',' '
' | grep -nx 'z' | cut -d: -f1)
    if [ -z "$z_column" ]; then
        echo "Error: No 'z' column found in $tx_file."
        exit 1
    fi
    z_values=$(zcat "$tx_file" | awk -F',' -v col="$z_column" 'NR>1 {print $col}' | sort -u)
    n_layer=$(echo "$z_values" | wc -l)
    echo "Number of voxel layers determined: $n_layer"
fi

# Run proseg with the specified parameters
echo "Running proseg with the specified parameters..."
proseg --cosmx \
    --nthreads $n_threads \
    --x-column x_global_px --y-column y_global_px \
    --output-expected-counts "$dir_output/cell-expected-counts.csv.gz" \
    --output-cell-metadata "$dir_output/cell-metadata.csv.gz" \
    --output-transcript-metadata "$dir_output/transcript-metadata.csv.gz" \
    --output-gene-metadata "$dir_output/gene-metadata.csv.gz" \
    --output-cell-polygons "$dir_output/cell-polygons.geojson.gz" \
    --output-cell-polygon-layers "$dir_output/cell-polygons-layers.geojson.gz" \
    --output-cell-hulls "$dir_output/cell-hulls.geojson.gz" \
    --output-cell-voxels "$dir_output/cell-voxels.csv.gz" \
    --voxel-layers "$n_layer" \
    "$tx_file"

# Notify user of successful completion
echo "Proseg processing completed successfully. Output files are stored in $dir_output."