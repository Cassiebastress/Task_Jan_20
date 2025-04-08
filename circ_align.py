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

    matches = 0
    for i in range(length):
        if seq1[i] != '-' and seq2[i] != '-' and seq1[i] == seq2[i]:
            matches += 1
    score = matches / length

    return round(score, 2)
