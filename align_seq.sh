#!/bin/bash

# Should take as input the fasta file to be aligned and create an aligned version with MAFFT

# First make sure there is a first parameter and the input file has .fasta extension
if [[ -z "$1" || ! "$1" == *.fasta ]]; then
	echo "Usage: $0 <input_fasta_file>"
	exit 1
fi

#Extract input filename
input_file="$1"

# Set extension to clustal format if the user wants
echo "Output in clustal format instead of default fasta format? (y/n)"
read output_format

# Check the user's input
if [[ "$ouput_format" == "y" || "$ouput_format" == "Y" ]]; then
	# Remove ".fasta" and add "_aln.fasta"
	output_file="${input_file%.fasta}_aln.fasta"

	# Run the mafft command with the auto alignment method
	mafft --auto "$input_file" > "$output_file"
else
	output_file="${input_file%.fasta}.aln"
	
	# Run mafft command with clustal formatted output
	mafft --auto --clustalout "$input_file" > "$output_file"
fi


