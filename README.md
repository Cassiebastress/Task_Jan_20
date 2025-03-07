# Task_Jan_20
 A small exercise to get familiar with how sequences and alignments are stored / represented

## align_seq.sh Script Overview
* Align sequences by passing the a fasta file to be aligned as the first argument
* The program will ask the user to choose between fasta and clustal format
* A file is created for the aligned sequence in their format of choice
* The file will be have the same base name as their input file
 * _aln.fasta will be appended if they requested fasta format
 * .aln will be appended if they requested clustal format

## Demo Files
* Three sequence files appear in the repository:
 1. FOXP2.fasta
    A file of sequences downloaded from the NCBI
 3. FOXP2_aln.fasta
    The output of running "align.seq.sh FOXP2" and not opting for clustal format
 5. FOXP2.aln
    The output of running "align.seq.sh FOXP2" and opting for clustal format

## Storage of Fasta and clustal files
1. Fasta
   * Starts with a header line which always begins with > followed by an identifier
   * Then the sequence
2. Clustal
   * Header line which describes the alignment
   * An identifier at the beginning of each line, and the aligned sequences on the right
   * Each column represents the same position in all of the sequences
   * The --- represents missing data in those positons of a sequence
   * The goal is to identify conserved regions and variable regions

## Setup Instructions
### 1. Python environment setup and installing dependencies
To keep dependencies isolated, create a virtual environment:

```bash
# Create a virtual environment in the .venv folder
python -m venv .venv

# activate the virtual environment
source .venv/bin/activate

# Install a package when you are in the virtual env
pip install -r requirements.txt

# deactivate it
deactivate
```

### 2. MAFFT Installation
MAFFT is required to use the sequence alignment tools in this project.
You can install MAFFT on macOS via Homebrew:

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install MAFFT
brew install mafft
```





