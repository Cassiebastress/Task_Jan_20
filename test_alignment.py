# Imports
import unittest
import os
from issue3 import align_fasta_to_seqs
from unittest.mock import patch
import subprocess
import tempfile
from Bio.Seq import Seq


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

    def test_mafft_not_installed(self):
        # Test that the function raises a RuntimeError
        # when MAFFT is not installed
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, "mafft")
            with self.assertRaises(RuntimeError) as e:
                align_fasta_to_seqs("input.fasta")
            self.assertIn("MAFFT is not installed", str(e.exception))

    def test_valid_input_file(self):
        # Create a temporary file with valid sequences
        with tempfile.NamedTemporaryFile("w", delete=False) as temp_file:
            temp_file.write(">TestSeq1\nACGTACGT\n>TestSeq2\nACGTACGA\n")
            temp_file_name = temp_file.name
        # Call the function with the temporary file
        try:
            result = align_fasta_to_seqs(temp_file_name)
            # Check that the result is a list
            self.assertIsInstance(result, list)
            # Check that list is not empty
            self.assertGreater(len(result), 0)
            # Check that the list contains Seq objects
            for seq in result:
                self.assertIsInstance(seq, Seq)
        # Remove the temporary file
        finally:
            os.remove(temp_file_name)


# Used to run the test suite when the scipt is executed directly
if __name__ == '__main__':
    unittest.main()
