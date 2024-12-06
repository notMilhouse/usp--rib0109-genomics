#!/bin/bash

# Set up folder structure to receive data
mkdir -p data/clinical
mkdir -p data/maf/jovem
mkdir -p data/maf/naojovem
mkdir -p data/rnaseq/jovem
mkdir -p data/rnaseq/naojovem

# Input file containing mappings (manifest and path)
INPUT_FILE="manifest_map"

# Command to run
COMMAND="gdc-client download"

# Check if the input file exists
if [[ ! -f $INPUT_FILE ]]; then
    echo "Error: File '$INPUT_FILE' not found."
    exit 1
fi
echo "downloading data described in manifests"
# Loop through each line of the file
while IFS= read -r line; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Extract the manifest and path from the line
    manifest=$(echo "$line" | awk '{print $1}')
    path=$(echo "$line" | awk '{print $2}')
    echo "Downloading data from '$manifest'"
    # Check if both values are present
    if [[ -z "$manifest" || -z "$path" ]]; then
        echo "Warning: Skipping invalid line '$line'"
        continue
    fi

    # Run the command with the manifest and path
    $COMMAND -m "$manifest" -d "$path"
done < "$INPUT_FILE"

echo "Finished processing all mappings."
