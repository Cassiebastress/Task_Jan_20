from issue3 import align_seqs
from Bio.Seq import Seq
from typing import List


def permutateSeq(seq: Seq, places: int) -> Seq:
    """
    Permutate a sequence by a given number of places
    """
    seq = list(seq)
    newSeq = []
    for i in range(len(seq)):
        newSeq.append(seq[(i + places) % len(seq)])
    return Seq("".join(newSeq))


def score_alignment(seq1: str, seq2: str) -> float:
    """
    Score an alignment
    """
    if (len(seq1) != len(seq2)):
        raise ValueError("Sequences must be the same length")
    if (len(seq1) == 0):
        raise ValueError("Sequences cannot be empty")
    length = len(seq1)
    seq1 = list(seq1)
    seq2 = list(seq2)

    matches = 0
    for i in range(length):
        if seq1[i] != '-' and seq2[i] != '-' and seq1[i] == seq2[i]:
            matches += 1
    score = matches / length

    return round(score, 2)


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
        cur_permutation = permutateSeq(seq1, i)
        aligned_seqs = align_seqs([cur_permutation, seq2])
        score = score_alignment(aligned_seqs[0], aligned_seqs[1])
        if score == 1.0:
            return aligned_seqs
        elif score > highest_score:
            highest_score = score
            best_alignment = aligned_seqs
    return best_alignment
