from issue3 import align_seqs


def permutateSeq(seq, places):
    """
    Permutate a sequence by a given number of places
    """
    seq = list(seq)
    newSeq = []
    for i in range(len(seq)):
        newSeq.append(seq[(i + places) % len(seq)])
    return "".join(newSeq)


def score_alignment(seq1, seq2):
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
    count = 0
    for i in range(length):
        if seq1[i] == "-" or seq2[i] == "-":
            count += 1
    score = (length - count) / length
    return round(score, 2)
