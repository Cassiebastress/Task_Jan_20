def still_has_content(data: list[list[str]], i: int) -> bool:
    """Return True if any template string has length > i."""
    return any(len(row[0]) > i for row in data)


def has_gap_at_position(data: list[list[str]], i: int) -> bool:
    """Return True if any alignment template has a '-' at position i."""
    return any(len(row[0]) > i and row[0][i] == '-' for row in data)


def find_sequence_with_char_at_i(data: list[list[str]], i: int) -> str:
    """Return the first character found at position i in any template."""
    for row in data:
        if len(row[0]) > i:
            return row[0][i]


def aligned_tuples_to_MSA(input_list: list[list[str]]) -> list[list[str]]:
    """
    Convert a list of single alignments into a multisequence alignment.
    """
    i = 0
    msa_ref = ""  # The final reference alignment string
    nonref = ['' for _ in input_list]  # Accumulate the aligned sequences

    while still_has_content(input_list, i):

        # If there's a gap in any template at position i
        if has_gap_at_position(input_list, i):
            # Add a gap to the reference sequence
            msa_ref += '-'

            for j in range(len(input_list)):
                template = input_list[j][0]
                sequence = input_list[j][1]

                # If the template has a gap at position i,
                # add the corresponding character from the sequence
                # to the nonref sequence
                if len(template) > i and template[i] == '-':
                    nonref[j] += sequence[i]

                # If the template has no gap at position i,
                # add a gap to the nonref sequence
                else:
                    nonref[j] += '-'
                    # Insert '-' in both template and sequence at position i
                    # to maintain alignment with other templates
                    input_list[j][0] = template[:i] + '-' + template[i:]
                    input_list[j][1] = sequence[:i] + '-' + sequence[i:]

        # If there's no gap in any template at position i
        else:
            # Take the character from any template to add to the reference
            char = find_sequence_with_char_at_i(input_list, i)
            msa_ref += char

            for j in range(len(input_list)):
                template = input_list[j][0]
                sequence = input_list[j][1]

                # If the template has a character at position i,
                # add the corresponding character from the sequence
                # to the nonref sequence
                if len(template) > i:
                    nonref[j] += sequence[i]

                # If the template has no character at position i,
                # add a gap to the nonref sequence
                # and pad the template with a gap to maintain alignment
                else:
                    nonref[j] += '-'
                    template += '-'

        i += 1

    # Return reference and aligned sequences as a list of strings
    return [msa_ref] + nonref
