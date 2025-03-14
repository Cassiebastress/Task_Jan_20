# Imports
import unittest
import os
from issue3 import align_fasta_to_seqs

# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio import SeqIO
# import tempfile
# import subprocess
# from typing import List


class TestSequenceAlignment(unittest.TestCase):

    def test_missing_input_file(self):
        # Test that the function raises
        # a FileNotFoundError when the input file
        # does not exist
        with self.assertRaises(FileNotFoundError) as e:
            align_fasta_to_seqs("non_existent_file.fasta")
        # Check that the error message contains the expected text
        self.assertIn("does not exist", str(e.exception))

    def test_invalid_character_in_input_file(self):
        # Create a file with invalid characters
        invalid_file = "invalid.fasta"
        try:
            with open(invalid_file, "w") as f:
                f.write(">TestSeq\nACGTNOTGOODCHARACTER\n")
            # Test that the function raises an error
            with self.assertRaises(ValueError) as e:
                # Call the function with the invalid file
                align_fasta_to_seqs(invalid_file)
            # Check that the error message contains the expected text
            self.assertIn("invalid characters", str(e.exception))
        # Remove the file
        finally:
            os.remove(invalid_file)


# Used to run the test suite when the scipt is executed directly
if __name__ == '__main__':
    unittest.main()
