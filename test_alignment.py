# Imports
import unittest
import os
from issue3 import align_fasta_to_seqs, align_seqs
from unittest.mock import patch
import subprocess
import tempfile
from Bio.Seq import Seq


class TestSequenceAlignment(unittest.TestCase):
    def test_missing_input_file(self):
        """
        Test that align_fasta_to_seqs raises a FileNotFoundError
        when the input file does not exist
        """
        with self.assertRaises(FileNotFoundError) as e:
            align_fasta_to_seqs("non_existent_file.fasta")
        self.assertIn("does not exist", str(e.exception))

    def test_invalid_character_in_input_file(self):
        """
        Test that align_fasta_to_seqs raises a ValueError
        when the input file contains invalid characters
        """
        invalid_file = "invalid.fasta"
        try:
            with open(invalid_file, "w") as f:
                f.write(">TestSeq\nACGTNOTGOODCHARACTER\n")
            with self.assertRaises(ValueError) as e:
                align_fasta_to_seqs(invalid_file)
            self.assertIn("invalid characters", str(e.exception))
        finally:
            os.remove(invalid_file)

    def test_mafft_not_installed(self):
        """
        Test that align_fasta_to_seqs raises a RuntimeError
        when MAFFT is not installed
        """
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, "mafft")
            with self.assertRaises(RuntimeError) as e:
                align_fasta_to_seqs("input.fasta")
            self.assertIn("MAFFT is not installed", str(e.exception))

    def test_valid_input_file(self):
        """
        Test that align_fasta_to_seqs returns a list of Seq objects
        when the input file is valid
        """
        with tempfile.NamedTemporaryFile("w+") as temp_file:
            temp_file.write(">TestSeq1\nACGTACGT\n>TestSeq2\nACGTACGA\n")
            temp_file.flush()  # Ensure the data is written before reading
            temp_file_name = temp_file.name
            result = align_fasta_to_seqs(temp_file_name)
            self.assertIsInstance(result, list)
            self.assertGreater(len(result), 0)
            for seq in result:
                self.assertIsInstance(seq, Seq)

    def test_empty_input_file(self):
        """
        Test that align_fasta_to_seqs raises a RuntimeError
        when the input file is empty
        """
        with tempfile.NamedTemporaryFile("w+") as temp_file:
            temp_file_name = temp_file.name
            with self.assertRaises(RuntimeError) as e:
                align_fasta_to_seqs(temp_file_name)
            self.assertIn("did not produce any output", str(e.exception))

    def test_single_valid_sequence(self):
        """
        Test that align_fasta_to_seqs returns a list of length 1
        when the input file contains a single valid sequence
        """
        with tempfile.NamedTemporaryFile("w+") as temp_file:
            temp_file.write(">TestSeq\nACGTACGT\n")
            temp_file.flush()
            result = align_fasta_to_seqs(temp_file.name)
            self.assertEqual(len(result), 1)
            self.assertIsInstance(result[0], Seq)

    def test_align_seqs_with_valid_input(self):
        seq_list = [Seq("ACGTACGT"), Seq("ACGTACGA")]
        result = align_seqs(seq_list)
        self.assertEqual(len(result), 2)
        for seq in result:
            self.assertIsInstance(seq, Seq)


# Used to run the test suite when the scipt is executed directly
if __name__ == '__main__':
    unittest.main()
