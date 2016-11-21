import sys
import argparse
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scipy.spatial import distance

from utilities.ev_couplings_v4 import EVcouplings
from ev_couplings_normalized import NormalizedEVcouplings


class SequenceProfile(object):
    """Representation of a coevolution-based profile of a specific sequence.

    Attributes:
        seq_id -- id of the sequence
        seq -- sequence
        mapped_seq -- of length N, seq mapped to target seq of the ev model
        seq_eij -- NxN matrix of eij values specific to the sequence
        profile -- coevolution-based profile of length N
    """

    def __init__(self, seq_id, seq, mapped_seq=None, seq_eij=None, profile=None):
        self.seq_id = seq_id
        self.seq = seq
        self.mapped_seq = mapped_seq
        self.seq_eij = seq_eij
        self.profile = profile

    def __repr__(self):
        return 'SequenceProfile(seq_id=%r, seq=%r, profile=%r)'\
            % (self.seq_id, self.seq, self.profile)

    @classmethod
    def from_repr(cls, obj_repr):
        """Construct SequenceProfile object from string representation."""
        return eval(obj_repr)

    @classmethod
    def create(cls, seq_record, ev_couplings):
        """Calculate profile of a given sequence."""
        seq_profile = cls(seq_record.id, str(seq_record.seq))
        seq_profile.map_to_model(ev_couplings)
        seq_profile.extract_seq_specific_eijs(ev_couplings)
        seq_profile.calc_profile()
        return seq_profile

    def map_to_model(self, ev_couplings):
        """Map sequence of interest to target sequence of the ev model."""
        aln_target_seq, aln_seq, _, _, _ = pairwise2.align.globalds(
            ev_couplings.target_seq.tostring().upper(), self.seq.upper(),
            blosum62, -11, -1
        )[0]
        self.mapped_seq = ''.join([res for i, res in enumerate(aln_seq)
                                   if not aln_target_seq[i] == '-'])

    def extract_seq_specific_eijs(self, ev_couplings):
        """Create matrix of eij values regarding a specific sequence."""
        aa_map = ev_couplings.alphabet_map
        eijs = ev_couplings.normalized_e_ij\
            if isinstance(ev_couplings, NormalizedEVcouplings)\
            else ev_couplings.e_ij
        self.seq_eij = np.zeros((len(self.mapped_seq), len(self.mapped_seq)))
        for i in xrange(len(self.mapped_seq)):
            for j in xrange(i+1, len(self.mapped_seq)):
                res_i, res_j = self.mapped_seq[i], self.mapped_seq[j]
                self.seq_eij[i, j] = self.seq_eij[j, i]\
                    = eijs[i, j, aa_map[res_i], aa_map[res_j]]

    def calc_profile(self):
        """Calculate profile."""
        self.profile = list(np.sum(self.seq_eij, axis=1))

    def dist(self, other):
        """Returns euclidean distance to another profile."""
        return distance.euclidean(self.profile, other.profile)


class EVprofiles(object):
    """Pairwise comparison of a query sequence to a set of other sequences
    by coevolution-based profiles.

    Attributes:
        query -- SeqRecord of the query sequence
        query_profile -- profile of the query sequence
        human_profiles -- other profiles the query gets compared to
                          (not necessarily human of course)
        dists_to_query -- distances of the human_profiles to query_profile
                          as list of tuples (human profile, distance to query)
    """

    # path to precalculated profiles
    HUMAN_PROFILES_FN = 'data/human_domains_from_scop.profiles'

    def __init__(self, query, eij_filename, human_seqs=None,
                 normalization_mode=NormalizedEVcouplings.ABS_MAX):
        ev_couplings = NormalizedEVcouplings(eij_filename, mode=normalization_mode)\
            if normalization_mode else EVcouplings(eij_filename)
        self.query = query
        self.query_profile = SequenceProfile.create(self.query, ev_couplings)
        if human_seqs is None:
            self.human_profiles = EVprofiles.load_profiles_from_file()
        else:
            self.human_profiles = [
                SequenceProfile.create(seq_record, ev_couplings)
                for seq_record in human_seqs
            ]
        self.dists_to_query = [
            (human_profile, self.query_profile.dist(human_profile))
            for human_profile in self.human_profiles
        ]
        self.dists_to_query = sorted(self.dists_to_query, key=lambda x: x[1])

    def to_file(self, out, num=1):
        """
            Write first num sequences with minimum distances from the query
            to file in fasta format.
        """
        top_dists_to_query = self.dists_to_query[:num]
        with open(out, 'w') as f:
            for profile, dist in top_dists_to_query:
                f.write('>%s distance_from_%s=%f\n' % (profile.seq_id, self.query.id, dist))
                f.write(profile.seq + '\n')
        print('resulting sequences written to %s' % out)

    @staticmethod
    def load_profiles_from_file():
        """Load precalculated profiles of a set of sequences from file."""
        with open(EVprofiles.HUMAN_PROFILES_FN, 'rU') as in_file:
            human_profiles = [
                SequenceProfile.from_repr(seq_profile_repr)
                for seq_profile_repr in in_file.readlines()
            ]
        return human_profiles


def init_parameters(args):
    """Read in and process parameters obtained from the command line."""
    with open(args.query) as query_file:
        query = SeqIO.read(query_file, 'fasta')
    try:
        with open(args.human_sequences, 'rU') as seq_file:
            human_seqs = list(SeqIO.parse(seq_file, 'fasta'))
    except TypeError:
        human_seqs = None
    return query, args.eij_file, human_seqs


def command_line():
    parser = argparse.ArgumentParser(
        description='Compare query sequence to a set of human sequences ' +
                    'by coevolution-based profiles'
    )
    parser.add_argument('--query', '-q', required=True,
                        help='Query sequence (in fasta format)')
    parser.add_argument('--eij_file', '-e', required=True,
                        help='eij file')
    parser.add_argument('--out', '-o', required=True,
                        help='Output file')
    parser.add_argument('--human_sequences', '-s', required=False,
                        help='Human sequences in fasta format ' +
                             '(default: precalculated profiles of human sequences are used)')
    parser.add_argument('--num_of_resulting_sequences', '-n', required=False,
                        help='The given number of sequences with minimal distance ' +
                             'to query will be printed to fasta file (default: 1)')
    args = parser.parse_args()

    ev_profiles = EVprofiles(*init_parameters(args))

    if args.num_of_resulting_sequences is not None:
        ev_profiles.to_file(args.out, int(args.num_of_resulting_sequences))
    else:
        ev_profiles.to_file(args.out)

if __name__ == '__main__':
    command_line()
