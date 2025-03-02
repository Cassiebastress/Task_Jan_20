# Demonstrating use of Seq in Biopython.

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import os
import subprocess

# FUNCTIONS:

# Function 1:
# Take input file and output file names
# Use Mafft to align sequences


def mafft_align(input_fasta, output_fasta):
    with open(output_fasta, "w") as outfile:
        subprocess.run(
            ["mafft", "--auto", input_fasta],
            stdout=outfile, stderr=subprocess.DEVNULL)


# Function 2:
# Take input file name
# Use mafft to align sequences
# Return list of Seq objects


def file_to_seq_list(input_fasta):
    # Use tempfile to create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a temporary output file
        temp_output = os.path.join(temp_dir, "aligned.fasta")

        # Use Function 1 function to align the sequences
        mafft_align(input_fasta, temp_output)

        # Read the aligned sequences and return a list of Seq objects
        seq_list = []
        for seq_record in SeqIO.parse(temp_output, "fasta"):
            seq_list.append(seq_record.seq)

    return seq_list

# Function 3:
# Take a list of Seq objects
# Return a list of Seq containing the allignments


def seq_list_to_seq_list(seq_list):
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a temporary input file
        temp_input = os.path.join(temp_dir, "input.fasta")

        # Write the sequences to a temporary FASTA file
        records = []
        for i, seq in enumerate(seq_list):
            record = SeqRecord(seq)
            record.id = "seq" + str(i)
            record.description = ""
            records.append(record)

        # Write the records to the temporary input file
        with open(temp_input, "w") as infile:
            SeqIO.write(records, infile, "fasta")

        # Use Function 3 to align the sequences and return a Seq list
        return file_to_seq_list(temp_input)

# TESTING:

# Function 1:
mafft_align("input.fasta", "output.fasta")

# Function 2:
aligned_seq_list_1 = file_to_seq_list("input.fasta")
for seq in aligned_seq_list_1:
    print(seq)

# Function 3:
test_seqs = [
    Seq("ACGTACGTACGT"),
    Seq("ACGTTGCAACGT"),
    Seq("ACGTACG"),
    Seq("TGCAACGT"),
]
aligned_seq_list_2 = seq_list_to_seq_list(test_seqs)
for seq in aligned_seq_list_2:
    print(seq)
