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
    Returns the number of matches between two sequences
    and the number of gap sections total for the two sequences ("breaks")
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

    # Check for breaks in seq1
    breaks = 0
    in_gap = False
    for i in range(length):
        if seq1[i] == '-':
            if not in_gap:
                breaks += 1
                in_gap = True
        else:
            in_gap = False
    # Check for breaks in seq2
    in_gap = False
    for i in range(length):
        if seq2[i] == '-':
            if not in_gap:
                breaks += 1
                in_gap = True
        else:
            in_gap = False

    return matches, breaks


def circular_align(seq1: Seq, seq2: Seq) -> List[Seq]:
    """
    Return the best alignment of two sequences
    """
    # if (len(seq1) != len(seq2)):
    #     raise ValueError("Sequences must be the same length")
    if (len(seq1) == 0 | len(seq2) == 0):
        raise ValueError("Sequences cannot be empty")
    highest_score = 0
    fewest_breaks = float('inf')
    best_alignment = [seq1, seq2]  # default to the original sequences
    for i in range(len(seq1)):
        cur_permutation = permutate_seq(seq1, i)
        aligned_seqs = align_seqs([cur_permutation, seq2])
        matches, breaks = score_alignment(aligned_seqs[0], aligned_seqs[1])
        if matches == len(seq1) | matches == len(seq2):
            return aligned_seqs
        elif matches > highest_score:
            highest_score = matches
            best_alignment = aligned_seqs
        elif matches == highest_score:
            # Check for fewer breaks if matches are tied
            if breaks < fewest_breaks:
                fewest_breaks = breaks
                best_alignment = aligned_seqs
    return best_alignment
