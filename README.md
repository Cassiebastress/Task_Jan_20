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
