"""Tool box solving miscellaneous tasks."""

import sys
import contextlib


@contextlib.contextmanager
def smart_open(filename=None):
    """Open file handle to write to or write to stdout.

    Parameters
    ----------
    filename : str
        Name of the output file. If `None` or '-', stdout is used.

    Yields
    ------
    fh : file
        Opened file handle or stdout.

    """
    fh = open(filename, 'w') if filename and filename != '-' else sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def split_args(parsed_args, *args, **kwargs):
    """Split args parsed from the command line in two sub dictionaries.

    Parameters
    ----------
    args : `Namespace`
        Arguments parsed by `argparse.ArgumentParser` from the command line.
    *args
        Variable length argument list containing keys.
    **kwargs
        Arbitrary keyword arguments. Supported keyword='this': (1) If `None`,
        a tuple of both generated dictionaries is returned. (2) If True,
        only the dictionary with matching keys as provided in `*args` is
        returned. (3) If false, the other dictionary is returned.

    Returns
    -------
    {tuple of dict, dict}
        Either both generated dictionaries as tuple are returned or a single
        dictionary, dependend on the the value of `this`.

    """
    sub_args_this = {k: v for k, v in vars(parsed_args).iteritems() if k in args}
    sub_args_other = {k: v for k, v in vars(parsed_args).iteritems() if k not in args}

    this = kwargs.get('this', None)
    if this is None:
        return sub_args_this, sub_args_other
    elif this:
        return sub_args_this
    else:
        return sub_args_other
