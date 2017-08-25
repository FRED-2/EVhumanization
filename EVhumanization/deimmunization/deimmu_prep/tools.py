"""Tool box for preparing and generating input files used to de-immunize an
amino acid sequence.

"""


def to_data_format(key_word, identifier, content, end=True, new_line=True):
    """Generalization and construction of an entry in data file.

    Parameters
    ----------
    key_word : {'param', 'set'}
        Start word of a data format entry.
    identifier : str
        Identifier for the specfication (just before ':=').
    content : str
        Content of the entry (right after ':=').
    end : bool, optional (default: True)
        If true, ';' is printed at the end of the entry (after `content`).
        Else, ';' is omitted.
    new_line : bool, optional (default: True)
        If true, a new line is printed after the entry.

    Returns
    -------
    str
        Single entry of data file.

    """
    stop = ';' if end else ''
    new = '\n' if new_line else ''
    return '%s %s := %s%s%s' % (key_word, identifier, content, stop, new)
