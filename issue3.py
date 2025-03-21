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
    with open(input_fasta, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
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
        with open(temp_output, "r") as infile:
            for seq_record in SeqIO.parse(infile, "fasta"):
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
