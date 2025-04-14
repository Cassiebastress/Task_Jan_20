import unittest
from circ_align import circular_align, permutate_seq, score_alignment
from Bio.Seq import Seq


class TestCircularAlign(unittest.TestCase):
    def test_permutate_seq(self):
        """
        Test the permutate_seq function
        """
        seq = Seq("ACGT")
        places = 2
        expected = Seq("GTAC")
        result = permutate_seq(seq, places)
        self.assertEqual(result, expected)

    def test_score_alignment(self):
        """
        Test the score_alignment function for partially aligned sequences
        """
        seq1 = Seq("ACGT-GTA")
        seq2 = Seq("ACGTTGTA")
        expected_score = 7
        result = score_alignment(seq1, seq2)
        self.assertEqual(result, expected_score)

    def test_score_empty_sequence(self):
        """
        Test that score_alignment raises a ValueError when sequences are empty
        """
        seq1 = Seq("")
        seq2 = Seq("")
        with self.assertRaises(ValueError) as e:
            score_alignment(seq1, seq2)
        self.assertIn("Sequences cannot be empty", str(e.exception))

    def test_score_different_length_sequences(self):
        """
        Test that score_alignment raises a ValueError when sequences are of
        different lengths
        """
        seq1 = Seq("ACGT")
        seq2 = Seq("ACG")
        with self.assertRaises(ValueError) as e:
            score_alignment(seq1, seq2)
        self.assertIn("Sequences must be the same length", str(e.exception))

    def test_align_empty_sequence(self):
        """
        Test that circular_align raises a ValueError when sequences are empty
        """
        seq1 = Seq("")
        seq2 = Seq("")
        with self.assertRaises(ValueError) as e:
            circular_align(seq1, seq2)
        self.assertIn("Sequences cannot be empty", str(e.exception))

    def test_align_different_length_sequences(self):
        """
        Test that circular_align raises a ValueError when sequences are of
        different lengths
        """
        seq1 = Seq("ACGT")
        seq2 = Seq("ACG")
        with self.assertRaises(ValueError) as e:
            circular_align(seq1, seq2)
        self.assertIn("Sequences must be the same length", str(e.exception))

    def test_circ_align(self):
        """
        Test the circular alignment of two sequences
        """
        seq1 = Seq("ACGTAAATTAAT")
        seq2 = Seq("AAACGTAAATTA")
        expected = [Seq("atacgtaaatta"), Seq("aaacgtaaatta")]
        result = circular_align(seq1, seq2)
        self.assertEqual(result, expected)

    def test_circ_align_perfect_alignment(self):
        """
        Test the circular alignment of two sequences with perfect match
        """
        seq1 = Seq("acgt")
        seq2 = Seq("acgt")
        expected = [seq1, seq2]
        result = circular_align(seq1, seq2)
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
