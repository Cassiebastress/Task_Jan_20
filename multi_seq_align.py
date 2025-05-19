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


def still_has_content(data: list[list[list[str]]], i: int) -> bool:
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

    while still_has_content(mutable_copy, i):
        if has_gap_at_position(mutable_copy, i):
            msa_ref += '-'
            for j, row in enumerate(mutable_copy):
                if len(row[0]) > i and row[0][i] == '-':
                    nonref[j].append(row[1][i])
                else:
                    nonref[j].append('-')
                    row[0].insert(i, '-')
                    row[1].insert(i, '-')
        else:
            char = find_sequence_with_char_at_i(mutable_copy, i)
            msa_ref += char
            for j, row in enumerate(mutable_copy):
                if len(row[0]) > i:
                    nonref[j].append(row[1][i])
                else:
                    nonref[j].append('-')
                    row[0].append('-')

        i += 1
    return [['ref', msa_ref]] + [[str(j + 1), ''.join(seq)] for j,
                                 seq in enumerate(nonref)]
