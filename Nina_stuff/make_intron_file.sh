#!/bin/bash

### The goal is to remove any lines that contain the word exon or is a header. 
# Prompt the user for the directory path and search string
read -p "Enter the directory path: " directory_path

# Check if the specified directory exists
if [ ! -d "$directory_path" ]; then
  echo "Error: Directory not found."
  exit 1
fi

# Prompt the user for the output directory path
read -p "Enter the output directory path: " output_directory_path

# Check if the output directory exists, if not, create it
if [ ! -d "$output_directory_path" ]; then
  mkdir -p "$output_directory_path"
fi

# Loop through all files in the directory
for file in "$directory_path"/*.vcf; do
  if [ -f "$file" ]; then
    # Get the filename without the directory path and extension
    filename=$(basename "$file" .vcf)
    
    # Use grep to find the lines containing the search string
    matching_lines=$(grep -vE "#|EXON" "$file")
    
    # Create the full path for the output file
    output_file="$output_directory_path/${filename}_no_header_or_exon.txt"
    
    # Save the matching lines to the output file
    echo "$matching_lines" > "$output_file"
    
    echo "Matching lines in '$file' saved to '$output_file'."
  fi
done
