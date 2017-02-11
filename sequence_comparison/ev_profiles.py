"""Module providing the functionality to calculate and compare
coevolution-based sequence profiles.

The profile of a sequence s is a numerical vector of the same length as s,
i.e. each profile value corresponds to a position in s.

A profile value p_i of position i in sequence s can be calculated in two
different ways:
    (1) Simple summation of e_ij pair couplings:
    p_i = sum_j=1^n e_ij(\sigma_i^s, \sigma_j^s)
    (2) Additional incorporating of h_i fields:
    p_i = h_i(\sigma_i^s) + 0.5 * sum_j=1^n e_ij(\sigma_i^s, \sigma_j^s)

"""

import sys
from argparse import ArgumentParser
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scipy.spatial import distance

from ev_couplings_normalized import ModificationMode
from utilities.ev_couplings_v4 import EVcouplings
from ev_couplings_normalized import ModifiedEVcouplings


class SequenceProfile(object):

    """Representation of the coevolution-based profile of a specific sequence.

    Parameters
    ----------
    seq_id : str
        ID of `seq`.
    seq : str
        Sequence for which the profile is to be calculated.
    profile : list of float, optional (default: None)
        The coevolution-based profile of `seq`.

    Attributes
    ----------
    seq_id : str
        ID of `seq`.
    seq : str
        Sequence for which the profile is to be calculated.
    mapped_seq : str (of length L)
        `seq` mapped to the target sequence of the model.
    seq_eij : np.array
        Matrix of shape (L, L) containing the e_ij pair couplings specific to
        `mapped_seq`.
    profile : np.array (of length L)
        The coevolution-based profile of `mapped_seq`.

    """

    def __init__(self, seq_id, seq, profile=None):
        self.seq_id = seq_id
        self.seq = seq
        self.profile = profile

    def __repr__(self):
        return ('SequenceProfile(seq_id=%r, seq=%r, profile=%r)'
                % (self.seq_id, self.seq, list(self.profile)))

    @classmethod
    def from_repr(cls, obj_repr):
        """Construct a sequence profile from string representation.

        Parameters
        ----------
        obj_repr : str
            String representation of a `SequenceProfile` object as generated
            by the function `__repr__()`.

        Returns
        -------
        `SequenceProfile`
            The sequence profile associated with the string representation
            provided.

        """
        profile = eval(obj_repr)
        profile.profile = np.array(profile.profile)
        return profile

    @classmethod
    def create(cls, seq_record, ev_couplings, incorporate_fields=True, top=None):
        """Calculate the profile of a given sequence.

        Parameters
        ----------
        seq_record : `SeqRecord`
            Sequence for which the profile is to be calculated.
        ev_couplings : `EVcouplings`
            e_ij binaries read in.
        incorporate_fields : bool, optional (default: True)
            If true, h_i fields are included in the profile calculation.
            Else, profiles are constructed only from e_ij pair couplings.
        top : {float, int}, optional (default: None)
            Number of top profile values used. All other profile positions are
            set to 0. If the value given is less than 1, then `top` will be
            interpreted as a fraction of all profile values. If `None`, all
            profile positions are used.

        Returns
        -------
        seq_profile : `SequenceProfile`
            The resulting sequence profile.

        """
        seq_profile = cls(seq_record.id, str(seq_record.seq))
        seq_profile.map_to_model(ev_couplings)
        seq_profile.extract_seq_specific_eijs(ev_couplings)
        seq_profile.calc_profile(ev_couplings, incorporate_fields, top)
        return seq_profile

    def map_to_model(self, ev_couplings):
        """Map the sequence of interest to the target sequence of the model.

        Parameters
        ----------
        ev_couplings : `EVcouplings`
            e_ij binaries read in.

        """
        aln_target_seq, aln_seq, _, _, _ = pairwise2.align.globalds(
            ev_couplings.target_seq.tostring().upper(), self.seq.upper(),
            blosum62, -11, -1
        )[0]
        self.mapped_seq = ''.join([res for i, res in enumerate(aln_seq)
                                   if not aln_target_seq[i] == '-'])

    def extract_seq_specific_eijs(self, ev_couplings):
        """Create the matrix of e_ij values specific to the mapped sequence.

        Parameters
        ----------
        ev_couplings : `EVcouplings`
            e_ij binaries read in.

        """
        aa_map = ev_couplings.alphabet_map
        eijs = ev_couplings.modified_e_ij\
            if isinstance(ev_couplings, ModifiedEVcouplings)\
            else ev_couplings.e_ij
        self.seq_eij = np.zeros((len(self.mapped_seq), len(self.mapped_seq)))
        for i in xrange(len(self.mapped_seq)):
            for j in xrange(i + 1, len(self.mapped_seq)):
                res_i, res_j = self.mapped_seq[i], self.mapped_seq[j]
                self.seq_eij[i, j] = self.seq_eij[j, i]\
                    = eijs[i, j, aa_map[res_i], aa_map[res_j]]

    def calc_profile(self, ev_couplings, incorporate_fields=True, top=None):
        """Calculate profile of the mapped sequence.

        Parameters
        ----------
        ev_couplings : `EVcouplings`
            e_ij binaries read in.
        incorporate_fields : bool, optional (default: True)
            If true, h_i fields are included in the profile calculation.
            Else, profiles are constructed only from e_ij pair couplings.
        top : {float, int}, optional (default: None)
            Number of top profile values used. All other profile positions are
            set to 0. If the value given is less than 1, then `top` will be
            interpreted as a fraction of all profile values. If `None`, all
            profile positions are used.

        """
        self.profile = np.zeros(len(self.mapped_seq))
        if incorporate_fields:
            for i, aa_i in enumerate(self.mapped_seq):
                self.profile[i] = (
                    ev_couplings.h_i[i, ev_couplings.alphabet_map[aa_i]] +
                    0.5 * np.sum(self.seq_eij[i])
                )
        else:
            self.profile = np.sum(self.seq_eij, axis=1)
        if top is not None:
            if top < 1:
                top *= ev_couplings.L
            profile_indexed = zip(range(len(self.profile)), self.profile)
            profile_sorted = sorted(profile_indexed, key=lambda x: x[-1],
                                    reverse=True)
            top_inds = [i for i, _ in profile_sorted[: int(top)]]
            for i, p_i in profile_indexed:
                self.profile[i] = p_i if i in top_inds else 0

    def dist(self, other):
        """Returns euclidean distance to another profile."""
        return distance.euclidean(self.profile, other.profile)


class EVprofiles(object):

    """Pairwise comparison of a query sequence to a set of other sequences
    by coevolution-based profiles.

    Parameters
    ----------
    query : `SeqRecord`
        The query sequence which gets compared to a set of other sequences.
    eij_filename : str
        Name of the binary file containing the model.
    human_seqs : list of `SeqRecord`, optional (default: None)
        The list of sequences that `query` gets compared to.
        If `None`, then precalculated profiles are used.
    modification_mode : str (ModificationMode.{ABS, ABS_MAX, MAX}), optional (default: ModificationMode.ABS)
        Mode to modify the e_ij pair couplings.
    incorporate_fields : bool, optional (default: True)
        If true, h_i fields are included in the profile calculation.
        Else, profiles are constructed only from e_ij pair couplings.
    top : {float, int}, optional (default: None)
        Number of top profile values used. All other profile positions are
        set to 0. If the value given is less than 1, then `top` will be
        interpreted as a fraction of all profile values. If `None`, all
        profile positions are used.

    Attributes
    ----------
    query : `SeqRecord`
        The query sequence which gets compared to a set of other sequences.
    query_profile : `SequenceProfile`
        The profile of `query`.
    human_profiles : list of `SequenceProfile`
        The list of profiles that `query_profile` gets compared to.
    dists_to_query : list of tuple of length 2 (`SequenceProfile`, float)
        List of euclidean distances from `query_profile` to `human_profiles`,
        as list of tuples: (human profile, distance to query profile).

    """

    # path to precalculated profiles
    HUMAN_PROFILES_FN = 'data/human_domains_from_digit.profiles'

    def __init__(self, query, eij_filename, human_seqs=None,
                 modification_mode=ModificationMode.ABS,
                 incorporate_fields=True, top=None):
        ev_couplings = ModifiedEVcouplings(eij_filename, mode=modification_mode)\
            if modification_mode else EVcouplings(eij_filename)

        # create query profile
        self.query = query
        self.query_profile = SequenceProfile.create(self.query, ev_couplings,
                                                    incorporate_fields, top)

        # get human profiles
        if human_seqs is None:
            self.human_profiles = EVprofiles.load_profiles_from_file()
        else:
            self.human_profiles = [
                SequenceProfile.create(seq_record, ev_couplings,
                                       incorporate_fields, top)
                for seq_record in human_seqs
            ]

        # calculate pairwise distances from the query to the human sequences
        np_human_profiles = np.empty((len(self.human_profiles), len(self.query_profile.profile)))
        for i, human_profile in enumerate(self.human_profiles):
            np_human_profiles[i] = human_profile.profile
        np_human_profiles -= self.query_profile.profile
        np_dists_to_query = np.linalg.norm(np_human_profiles, axis=1)
        self.dists_to_query = [(profile, np_dists_to_query[i])
                               for i, profile in enumerate(self.human_profiles)]
        self.dists_to_query = sorted(self.dists_to_query, key=lambda x: x[1])

    def to_file(self, out, num=1):
        """Write top human sequence(s) to fasta file.

        The given number of profiles with minimum distance to the query are
        written to file.

        Parameters
        ----------
        out : str
            Name of the ouput file.
        num : int
            Number of sequences to be written to file. If `None`, all sequences
            are written to file (in order of increasing distance to the query).

        """
        top_dists_to_query = self.dists_to_query[:num]\
            if num is not None else self.dists_to_query
        with open(out, 'w') as f:
            for profile, dist in top_dists_to_query:
                f.write('>%s distance_from_%s=%f\n'
                        % (profile.seq_id, self.query.id, dist))
                f.write(profile.seq + '\n')

    def to_stdout(self, num=1):
        """Write top human sequence(s) to stdout in fasta format.

        The given number of profiles with minimum distance to the query are
        written to stdout.

        Parameters
        ----------
        num : int
            Number of sequences to be written to file. If `None`, all sequences
            are written to file (in order of increasing distance to the query).

        """
        top_dists_to_query = self.dists_to_query[:num]\
            if num is not None else self.dists_to_query
        for profile, dist in top_dists_to_query:
            print ('>%s distance_from_%s=%f\n%s'
                   % (profile.seq_id, self.query.id, dist, profile.seq))

    @staticmethod
    def load_profiles_from_file():
        """Load precalculated profiles of a set of sequences from file."""
        with open(EVprofiles.HUMAN_PROFILES_FN, 'rU') as in_file:
            return [
                SequenceProfile.from_repr(seq_profile_repr)
                for seq_profile_repr in in_file.readlines()
            ]


def init_params(args):
    """Initialize arguments parsed from the command line."""
    args.query = SeqIO.read(args.query, 'fasta')
    args.human_seqs = args.human_seqs if args.human_seqs is None\
        else list(SeqIO.parse(args.human_seqs, 'fasta'))
    args.modification_mode = ModificationMode.getMode(args.modification_mode)\
        if args.modification_mode.lower() != 'none' else None
    return args


def command_line():
    """Define parser to read arguments from the command line."""
    parser = ArgumentParser(description="""Pairwise comparison of a query
                            sequence to a set of other sequences by
                            coevolution-based profiles.""")
    parser.add_argument('query', help="""The query sequence which gets compared
                        to a set of other sequences. Usually, this is the non-human
                        antibody sequence to be humanized.""")
    parser.add_argument('eij_filename', help="""e_ij binary file.""")
    parser.add_argument('--modification_mode', '-m', required=False,
                        choices=['none', 'max', 'abs_max', 'abs'],
                        default='abs', help="""Modification mode, applied to
                        the e_ij pair couplings. Must be one of: 'none' (no
                        modification), 'max' (normalize each e_ij matrix by its
                        absolute maximum value), 'abs_max' (normalize each absolute
                        e_ij matrix by its maximum value) or 'abs' (set each e_ij
                        entry to its absolute value). The default is 'abs'.""")
    parser.add_argument('--top', '-t', required=False, type=float,
                        help="""Number of top profile values used. All other
                        profile positions are set to 0. If the value given is
                        less than 1, then it will be interpreted as a fraction
                        of all profile values. As default, all profile positions
                        are used.""")
    parser.add_argument('--human_seqs', '-s', required=False,
                        help="""The sequences that the query sequence gets
                        compared to. Usually, this is a set of human sequences.
                        As default, precalculated profiles are used.""")
    parser.add_argument('--num_of_resulting_sequences', '-n', required=False,
                        type=int, default=1, help=""""The value given defines the
                        number of top sequences that will be printed. The N top
                        sequences are the N first sequences whose profiles have
                        minimal distance to the queries profile (default: N=1).""")
    parser.add_argument('--simple_eij_summation', dest='incorporate_fields',
                        action='store_false', help="""Calculate profiles without
                        incorporating fields. A profile value will be calculated
                        as a simple summation of e_ij pair couplings specific
                        for the sequence of interest.""")
    parser.add_argument('--out', '-o', required=False, help='Output file')
    return init_params(parser.parse_args())


def main(args):
    """Initialize EV profiles and write results."""
    ev_profiles = EVprofiles(args.query, args.eij_filename, args.human_seqs,
                             args.modification_mode, args.incorporate_fields,
                             args.top)
    ev_profiles.to_file(args.out, num=args.num_of_resulting_sequences)\
        if args.out is not None\
        else ev_profiles.to_stdout(num=args.num_of_resulting_sequences)


if __name__ == '__main__':
    args = command_line()
    main(args)
