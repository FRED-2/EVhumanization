import sys
from argparse import ArgumentParser
from Bio import SeqIO
import time

from ev_profiles import SequenceProfile
from ev_profiles import EVprofiles

from utilities.ev_couplings_v4 import EVcouplings
from ev_couplings_normalized import ModifiedEVcouplings, ModificationMode


def precalculate_profiles(human_seqs, eij_filename, out,
                          modification_mode=ModificationMode.ABS,
                          incorporate_fields=True):
    """Precalculate profiles and write their string representation to file.

    Parameters
    ----------
    human_seqs : list of `SeqRecord`
        Sequences for which profiles are to be precalculated.
    eij_filename : str
        Name of the binary e_ij file.
    out : str
        Name of the output file.
    mode : str (ModificationMode.{MAX, ABS_MAX, ABS}), optional (default: ModificationMode.ABS)
        Modification mode, applied to the e_ij pair couplings.
    incorporate_fields : bool, optional (default: True)
        If true, h_i fields are included in the profile calculation.
        Else, profiles are constructed only from e_ij pair couplings.

    """
    ev_couplings = EVcouplings(eij_filename) if modification_mode is None\
        else ModifiedEVcouplings(eij_filename, modification_mode)
    start = time.time()
    with open(out, 'w') as f:
        for i, seq_record in enumerate(human_seqs):
            f.write(repr(SequenceProfile.create(seq_record, ev_couplings)) + '\n')
            print >> sys.stderr, i + 1, '/', len(human_seqs), 'done'
    end = time.time()
    print >> sys.stderr, 'Elapsed time: %.2fs' % (end - start)


def init_params(args):
    """Initialize arguments parsed from the command line."""
    args.human_seqs = list(SeqIO.parse(args.human_seqs, 'fasta'))
    args.modification_mode = ModificationMode.getMode(args.modification_mode)\
        if args.modification_mode.lower() != 'none' else None
    return args


def command_line():
    """Define parser to read arguments from the command line."""
    parser = ArgumentParser(description="""Precalculate coevolution-based
                            sequence profiles of a set of sequences.""")
    parser.add_argument('human_seqs', help="""The sequences that the query
                        sequence gets compared to. Usually, this is a set of
                        human sequences. As default, precalculated profiles are
                        used.""")
    parser.add_argument('eij_filename', help="""e_ij binary file.""")
    parser.add_argument('out', help='Output file')
    parser.add_argument('--modification_mode', '-m', required=False,
                        choices=['none', 'max', 'abs_max', 'abs'],
                        default='abs', help="""Modification mode, applied to
                        the e_ij pair couplings. Must be one of: 'none' (no
                        modification), 'max' (normalize each e_ij matrix by its
                        absolute maximum value), 'abs_max' (normalize each absolute
                        e_ij matrix by its maximum value) or 'abs' (set each e_ij
                        entry to its absolute value). The default is 'abs'.""")
    parser.add_argument('--simple_eij_summation', dest='incorporate_fields',
                        action='store_false', help="""Calculate profiles without
                        incorporating fields. A profile value will be calculated
                        as a simple summation of e_ij pair couplings specific
                        for the sequence of interest.""")
    return init_params(parser.parse_args())


def main(args):
    precalculate_profiles(**vars(args))


if __name__ == '__main__':
    args = command_line()
    main(args)
