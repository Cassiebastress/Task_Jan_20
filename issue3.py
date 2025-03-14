# Demonstrating use of Seq in Biopython.

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import os
import subprocess
from typing import List

# FUNCTIONS:

# Function 1:
# Take input file and output file names
# Use Mafft to align sequences


def align_fasta_file(input_fasta: str, output_fasta: str) -> None:
    # Check that proper input_fasta was provided
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(
            f"Error: The input file '{input_fasta}' does not exist "
            "or is not a file."
        )

    # Check that input_fasta contains only A, C, G, T characters
    valid_chars = set("ACGT")
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_upper = set(record.seq.upper())
        if not seq_upper.issubset(valid_chars):
            raise ValueError(
                f"Error: Sequence '{record.id}' in the input file "
                f"'{input_fasta}' contains invalid characters."
            )

    # Check that MAFFT is installed and available in the Path
    try:
        subprocess.run(
            ["mafft", "--version"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError:
        raise RuntimeError(
            "Error: MAFFT is not installed or not available in the PATH."
        )

    # Run MAFFT and capture stderr
    with open(output_fasta, "w") as outfile:
        result = subprocess.run(
            ["mafft", "--auto", input_fasta],
            stdout=outfile, stderr=subprocess.PIPE, text=True)

        # Check if MAFFT failed:
        if result.returncode != 0:
            os.remove(output_fasta)
            raise RuntimeError(
                f"Error: MAFFT failed with error:\n {result.stderr}"
            )
        # Check if MAFFT produced empty output
        elif os.path.exists(output_fasta) and os.path.getsize(output_fasta) == 0:
            os.remove(output_fasta)
            raise RuntimeError(
                f"Error: MAFFT did not produce any output for '{input_fasta}'."
            )


# Function 2:
# Take input file name
# Use mafft to align sequences
# Return list of Seq objects


def align_fasta_to_seqs(input_fasta: str) -> List[Seq]:
    # Use tempfile to create a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a temporary output file
        temp_output = os.path.join(temp_dir, "aligned.fasta")

        # Use Function 1 function to align the sequences
        align_fasta_file(input_fasta, temp_output)

        # Read the aligned sequences and return a list of Seq objects
        seq_list = []
        for seq_record in SeqIO.parse(temp_output, "fasta"):
            seq_list.append(seq_record.seq)

    return seq_list

# Function 3:
# Take a list of Seq objects
# Return a list of Seq containing the allignments


def align_seqs(seq_list: List[Seq]) -> List[Seq]:
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
        return align_fasta_to_seqs(temp_input)

# # TESTING:

# # Function 1:
# align_fasta_file("input.fasta", "output.fasta")

# # Function 2:
# aligned_seq_list_1 = align_fasta_to_seqs("input.fasta")
# for seq in aligned_seq_list_1:
#     print(seq)


# # Function 3:
# test_seqs = [
#     Seq("ACGTACGTACGT"),
#     Seq("ACGTTGCAACGT"),
#     Seq("ACGTACG"),
#     Seq("TGCAACGT"),
# ]
# aligned_seq_list_2 = align_seqs(test_seqs)
# for seq in aligned_seq_list_2:
#     print(seq)

# # Improper Input Tests:

# # Function 1:
# align_fasta_file("not_a_file.fasta", "output.fasta")
# align_fasta_file("bad_input.fasta", "bad_output.fasta")  # Not a valid FASTA file

# # Function 2:
# align_fasta_to_seqs("not_a_file.fasta")
# align_fasta_to_seqs("bad_input.fasta")  # Not a valid FASTA file
