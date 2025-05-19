import copy


def prepare_mutable_copy(data: list[list[str]]) -> list[list[list[str]]]:
    """
    Create a deep copy of the data where the first string in each row is
    converted to a list of characters for mutability.

    Returns:
        A new list where both rows are lists of characters (mutable).
    """
    new_data = copy.deepcopy(data)
    for row in new_data:
        row[0] = list(row[0])
        row[1] = list(row[1])
    return new_data


def still_has_real_content(data: list[list[list[str]]], i: int) -> bool:
    """
    Returns True if any sequence has a non-dash character at or after index i.
    """
    for row in data:
        for j in range(i, len(row[0])):
            if len(row[0]) > j:
                return True
            # if row[0][j] != '-':
            #     return True
    return False


def has_gap_at_position(data: list[list[str]], i: int) -> bool:
    """
    Check if any of the first elements (index 0) in the inner lists
    has a '-' character at the given index i.

    Args:
        data (list): The 2D list of strings.
        i (int): The position to check for a gap ('-').

    Returns:
        bool: True if any first element has '-' at position i, False otherwise.
    """
    for row in data:
        if len(row[0]) > i and row[0][i] == '-':
            return True
    return False


def find_sequence_with_char_at_i(data: list[list[str]], i: int) -> str:
    """
    Find the first character at index i in the first elements of the inner
    lists.

    Args:
        data (list): The 2D list of strings.
        i (int): The index to check.

    Returns:
        str: The character found at index i, or '-' if not found.
    """
    for row in data:
        if len(row[0]) > i:
            return row[0][i]
    return '-'


def aligned_tuples_to_MSA(input_list: list) -> list:
    i = 0
    msa_ref = ""
    mutable_copy = prepare_mutable_copy(input_list)

    # Initialize nonref as a list of empty lists
    nonref = [[] for _ in range(len(mutable_copy))]

    while still_has_real_content(mutable_copy, i):
        if has_gap_at_position(mutable_copy, i):
            print(f"Gap found at position {i}")
            print("adding gap to reference")
            msa_ref += '-'
            for j, row in enumerate(mutable_copy):
                if len(row[0]) > i and row[0][i] == '-':
                    print(f"Ref sequence {j+1} has gap at position {i}")
                    print(f"Adding char {row[1][i]} to output sequence {j+1}")
                    nonref[j].append(row[1][i])
                else:
                    print(f"Sequence {j+1} has no gap at position {i}")
                    print(f"Adding '-' to output sequence {j+1}")
                    nonref[j].append('-')
                    print(f"Adding '-' to ref sequence {j+1} at position {i}")
                    row[0].insert(i, '-')
                    print(f"Adding '-' to nonref sequence {j+1} at position {i}")
                    row[1].insert(i, '-')
        else:
            print(f"No gap at position {i}")
            char = find_sequence_with_char_at_i(mutable_copy, i)
            print(f"Adding character '{char}' to reference")
            msa_ref += char
            for j, row in enumerate(mutable_copy):
                if len(row[0]) > i:
                    print(f"Nonref sequence {j+1} has character '{row[1][i]}' at position {i}")
                    print(f"Adding character '{row[1][i]}' to sequence {j+1}")
                    nonref[j].append(row[1][i])
                else:
                    print(f"Nonref sequence {j+1} has no character at position {i}")
                    print(f"Adding '-' to output sequence {j+1}")
                    nonref[j].append('-')
                    print(f"Adding '-' to reference sequence {j+1}")
                    row[0].append('-')

        i += 1
    return [['ref', msa_ref]] + [[str(j + 1), ''.join(seq)] for j,
                                 seq in enumerate(nonref)]
