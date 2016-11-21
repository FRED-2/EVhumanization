import sys
import argparse
from Bio import SeqIO

from ev_profiles import SequenceProfile
from ev_profiles import EVprofiles

from ev_couplings_v4 import EVcouplings
from utilities.ev_couplings_normalized import NormalizedEVcouplings


def precalculate_profiles(human_seqs, eij_file, mode=NormalizedEVcouplings.ABS_MAX):
    """Calculate profiles based on normalized eijs."""
    ev_couplings = EVcouplings(eij_file) if mode == 'no'\
        else NormalizedEVcouplings(eij_file, mode)
    return [SequenceProfile.create(seq_record, ev_couplings)
            for seq_record in human_seqs]


def write_to_file(profiles, out):
    """Write profiles to file."""
    with open(out, 'w') as out_file:
        for p in profiles:
            out_file.write(repr(p) + '\n')
    print 'profiles written to %s' % out


def command_line():
    parser = argparse.ArgumentParser(
        description='Precalculate coevolution-based profiles of given human sequences'
    )
    parser.add_argument('--sequences', '-s', required=True,
                        help='Sequences whose profiles are to be precalculated ' +
                             '(in fasta format)')
    parser.add_argument('--eij_file', '-e', required=True, help='eij file')
    parser.add_argument('--out', '-o', required=True, help='Output file')
    parser.add_argument('--normalization', '-n', required=False,
                        help="Normalization mode (one of: 'max', 'abs max', 'no'); " +
                             'no: no normalization, ' +
                             'max: normalize each eij matrix by its maximum value, ' +
                             'abs max: normalize each absolute eij matrix by its maximum value (default)')
    args = parser.parse_args()

    with open(args.sequences, 'rU') as seq_file:
        human_seqs = list(SeqIO.parse(seq_file, 'fasta'))

    profiles = precalculate_profiles(human_seqs, args.eij_file)\
        if args.normalization is None\
        else precalculate_profiles(human_seqs, args.eij_file, args.normalization)

    write_to_file(profiles, args.out)

if __name__ == '__main__':
    command_line()
