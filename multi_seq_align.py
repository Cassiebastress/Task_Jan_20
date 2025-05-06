def shift(char_list):
    """
    Shift the list of characters to the left by one position, in-place.
    """
    print(f"Before shift: {char_list}")
    if char_list:
        del char_list[0]
    print(f"After shift: {char_list}")


def find_true_length(input: str) -> int:
    """
    Find the true length of a sequence by removing gaps.

    Args:
        input (str): The input sequence.

    Returns:
        int: The true length of the sequence.
    """
    length = len([c for c in input if c != '-'])
    print(f"True length of sequence: {length}")
    return length


def shiftAll(input: list) -> None:
    """
    Shift all sequences in the input list to the left by one position.

    Args:
        input (list): The input list of sequences.

    Returns:
        list: The shifted list of sequences.
    """
    print(f"Shifting all sequences: {input}")
    print(f"length of input is: {len(input)}")
    for i in range(len(input)):
        shift(input[i][0])


def check_for_gaps(input: list) -> bool:
    print(f"Checking for gaps in input: {input}")
    for j in range(len(input)):
        # If the reference sequence has a gap
        if input[j][0][0] == '-':
            shift(input[j][0])
            return True
    return False


def aligned_tuples_to_MSA2(input: list) -> str:

    # Initialize variables
    ref = ""
    num = find_true_length(input[0][0])  # Use first reference sequence
    hasGap = False

    for i in range(num):
        hasGap = check_for_gaps(input)
        if hasGap:
            # while there is a gap in one of the ref seqs
            # add a gap to the final reference sequence
            # and shift that reference input seq to the left
            while (hasGap):
                print("Adding gap to reference")
                ref += '-'
                hasGap = check_for_gaps(input)
            # No longer any gaps, so add the first character of any
            # sequence to the reference sequence since all sequences
            # are aligned at this position
            print(f"Adding character to reference: {input[0][0][0]}")
            ref += input[0][0][0]
            shiftAll(input)
        else:
            # Add the first character of any sequence to the
            # reference sequence since all sequences are now aligned
            # at this position
            print(f"Adding character to reference: {input[0][0][0]}")
            ref += input[0][0][0]
            shiftAll(input)

    return ref
