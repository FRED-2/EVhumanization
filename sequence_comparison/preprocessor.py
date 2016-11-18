import sys
import argparse
from Bio import SeqIO

from ev_profiles import SequenceProfile
from ev_profiles import EVprofiles

sys.path.append('utilities')
from ev_couplings_v4 import EVcouplings


def init_parameters(args):
    """Read in parameters."""
    with open(args.sequences, 'rU') as seq_file:
        human_seqs = list(SeqIO.parse(seq_file, 'fasta'))
    return human_seqs, EVcouplings(args.eij_file)


def precalculate_profiles(human_seqs, ev_couplings):
    """Calculate profiles based on normalized eijs."""
    normed_eij = EVprofiles.normalize_eijs(ev_couplings.e_ij)
    return [SequenceProfile.create(seq_record, ev_couplings, normed_eij)
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
    args = parser.parse_args()

    profiles = precalculate_profiles(*init_parameters(args))
    write_to_file(profiles, args.out)

if __name__ == '__main__':
    command_line()
