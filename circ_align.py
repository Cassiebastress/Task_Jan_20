from issue3 import align_seqs
from Bio.Seq import Seq
from typing import List


def permutate_seq(seq: Seq, places: int) -> Seq:
    """
    Permutate a sequence by a given number of places
    """
    rotation = places % len(seq)
    return seq[-rotation:] + seq[:- rotation]


def score_alignment(seq1: str, seq2: str) -> float:
    """
    Score an alignment
    """
    if (len(seq1) != len(seq2)):
        raise ValueError("Sequences must be the same length")
    if (len(seq1) == 0):
        raise ValueError("Sequences cannot be empty")

    length = len(seq1)
    matches = 0
    for i in range(length):
        if seq1[i] != '-' and seq1[i] == seq2[i]:
            matches += 1

    return matches


def circular_align(seq1: Seq, seq2: Seq) -> List[Seq]:
    """
    Return the best alignment of two sequences
    """
    if (len(seq1) != len(seq2)):
        raise ValueError("Sequences must be the same length")
    if (len(seq1) == 0):
        raise ValueError("Sequences cannot be empty")
    highest_score = 0
    best_alignment = [seq1, seq2]  # default to the original sequences
    for i in range(len(seq1)):
        cur_permutation = permutate_seq(seq1, i)
        aligned_seqs = align_seqs([cur_permutation, seq2])
        score = score_alignment(aligned_seqs[0], aligned_seqs[1])
        if score == len(seq1):
            return aligned_seqs
        elif score > highest_score:
            highest_score = score
            best_alignment = aligned_seqs
    return best_alignment
