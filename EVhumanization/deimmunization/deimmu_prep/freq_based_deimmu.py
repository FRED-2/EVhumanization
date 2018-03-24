#!/usr/local/bin/python2.7
# encoding: utf-8

"""Preparation and generation of input files
used to de-immunize an antibody sequence.

Supported way of setting mutations:
Mutations to residues that are frequently seen
in that position in murine organisms.
"""
import os
import sys
import ConfigParser
import pickle

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))

from argparse import ArgumentParser

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from EVhumanization.deimmunization.grafting.kabat_numbering import KabatNumbering
from EVhumanization.utilities.ev_couplings_v4 import EVcouplings
from EVhumanization.utilities.smart_tools import split_args
from abstract_deimmu_prep import AbstractDeimmuPreparation


class AntibodyDeimmuPreparation(AbstractDeimmuPreparation):

    """Preparation and generation of input
    files used to de-immunize
    an antibody sequence.

    This class supports mutations to
    residues that are frequently seen
    in that position in murine organisms.

    Parameters
    ----------
    config : `ConfigParser`
        Config file, describing further parameters for de-immunization.
    alignment : `MultipleSeqAlignment`
        Alignment, used for plmc calculation in focus mode
        (with upper and lower letters).
    ev_couplings : `EVcouplings`
        e_ij binaries read in, generated from `alignment`.
    source_seq : str
        Sequence which serves as source of backmutation residues. Usually,
        this is the murine sequence to be humanized.
    include_similar_amino_acids : bool, optional (default: False)
        True, if amino acids similiar to the `source_seq` residues are
        also considered.

    Attributes
    ----------
    wildtype : `Wildtype`
        The wildtype to be de-immunized. This is usually a human sequence
        onto which murine residues (e.g. CDRs) are already grafted.
    num_mutations : int
        Number of mutations to be introduced in the wildtype.
    epitope_length : int
        Length of an epitope.
    allele_coll : `AlleleCollection`
        Collection of alleles.
    excluded_pos : list of int
        List of positions not allowed to be mutated.
    ignored_pos : list of int
        List of positions that are fully ignored.
    possible_mutations : list of set of str
        For each position, defines a set of amino acids allowed
        in that position.
    hi : dict (key=tuple of length 2 (int, str), value=float)
        h_i fields relevant to the calculation,
        as dictionary: key=(i, aa_i), value=h_i of amino acid aa_i in position i
    eij_indices : list of tuple of length 2 (int, int)
        List of paired e_ij indices whose associated values are relevant
        to the calculation
    eij : dict (key=tuple of length 4 (int, int, str, str), value=float)
        e_ij pair couplings relevant to the calculation,
        as dictionary: key=(i, j, aa_i, aa_j), value=e_ij of amino acids
        aa_i and aa_j in positions i and j, respectively.

    """

    def __init__(self, config, alignment, ev_couplings,
                 chain="heavy", freq_thresh=0.01):
        super(AntibodyDeimmuPreparation, self).__init__(config, alignment, ev_couplings)
        self.set_possible_mutations(chain=chain, freq_thresh=freq_thresh)
        self.extract_ev_paras(ev_couplings)

    def set_possible_mutations(self, chain="heavy", freq_thresh=0.01):
        """For each wildtype position, generate a set of amino acid
        substitutions (the wildtype amino acid included).

        Method fills the attribute `possible_mutations`.

        Parameters
        ----------
        chain : {"heavy", "light"}, optional (default: "heavy")
            Type of the variable antibody chain that is
            to be humanized. Can either be "heavy" or
            "light".
        freq_thresh : float, optional (default: 0.8)
            Amino acids in some position seen
            with a frequency >= the given threshold
            are considered as possible substitution
            in that position.

        """
        chain = chain.strip().lower()
        chain_key = chain[0].upper()

        # load murine amino acid frequencies
        # for every Kabat position from file
        filename = "EVhumanization/data/aa_freqs_{}_mus.pcl".format(chain)
        aa_freqs = pickle.load(open(filename, "rb"))

        # Kabat numbering of the sequence of interest
        kabat = KabatNumbering(
            SeqRecord(Seq(self.wildtype.sequence))
        ).kabat_dict

        kabat = {
            value: "".join(map(str, key))
            for key, value in kabat.iteritems()
        }

        print 'Possible mutations:'
        self.possible_mutations = [set(res) for res in self.wildtype.sequence]
        for i, muts in enumerate(self.possible_mutations):
            if i not in self.excluded_pos and i not in self.ignored_pos:
                pos = i + 1
                kabat_pos = kabat[pos]
                try:
                    freqs = aa_freqs[chain_key + kabat_pos]
                except KeyError as e:
                    print "No frequency information for position", e.message
                    continue
                above_thresh = filter(lambda t: t[-1] >= freq_thresh * 100, freqs)
                muts.update(map(lambda t: t[0], above_thresh))
            print muts


def init_params(args):
    """Initialize arguments parsed from the command line."""
    config = ConfigParser.ConfigParser()
    config.read(args.config)
    args.config = config

    args.alignment = AlignIO.read(args.alignment, 'fasta')
    ev_couplings = EVcouplings(args.eij_filename, file_format="plmc_v1")
    vars(args)['ev_couplings'] = ev_couplings
    del vars(args)['eij_filename']

    return args


def command_line():
    """Define parser to read arguments from the command line."""
    parser = ArgumentParser(description="""Preparation and generation of input
                            files used to de-immunize an amino acid sequence""")
    parser.add_argument('--config', '-c', required=True, help="""Config file,
                        describing further parameters for de-immunization.""")
    parser.add_argument('--eij_filename', '-e', required=True, help='e_ij binary file.')
    parser.add_argument('--alignment', '-a', required=True, help="""Alignment,
                        used for plmc calculation in focus mode (with upper and
                        lower letters). The sequence to be de-immunized must be
                        the first record in the alignment.""")
    parser.add_argument('--out', '-o', required=True, help="""The output file
                        of the data model (in AMPL/GMPL) format.""")
    parser.add_argument('--model', '-m', required=True, nargs=2, help="""The ILP
                        model files in MathProg (first argument: imm, second
                        argument: en).""")
    parser.add_argument("--chain", "-ch", required=True, choices=["heavy", "light"],
                        help=""""Type of variable antibody chain that is to be
                        humanized. Can either be 'heavy' or 'light'.""")
    parser.add_argument("--freq_thresh", "-f", default=0.01, type=float,
                        help=""""Amino acids with a frequency above this threshold
                        in a certain position will be considered as possible mutations.""")
    return init_params(parser.parse_args())


def main(args):
    """Call antibody de-immunization preparation and write data and lp files."""
    out_args, main_args = split_args(args, 'model', 'out')
    deimmu = AntibodyDeimmuPreparation(**main_args)
    deimmu.to_data_file(args.out)
    deimmu.generate_lp_files(**out_args)


if __name__ == '__main__':
    args = command_line()
    main(args)
