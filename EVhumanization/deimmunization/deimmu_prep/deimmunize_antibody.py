#!/usr/local/bin/python2.7
# encoding: utf-8
"""Preparation and generation of input files used to de-immunize an
antibody sequence.

Supported are two ways of setting possible mutations: (1) Only backmutating to
the corresponding murine amino acid is allowed. (2) Amino acids similiar to the
murine one are also considered.

"""
import os
import sys
import ConfigParser

sys.path.append(os.path.join(os.path.dirname(__file__), '../../..'))

from argparse import ArgumentParser

from Bio import AlignIO
from Bio import SeqIO

from EVhumanization.utilities.ev_couplings_v4 import EVcouplings
from EVhumanization.utilities.smart_tools import split_args
from abstract_deimmu_prep import AbstractDeimmuPreparation

_SMALL_HYDROPHOBIC_AA = list('AG')
_LARGE_HYDROPHOBIC_AA = list('VLIMPFW')
_POLAR_AA = list('STNQCY')
_POS_CHARGED_AA = list('KRH')
_NEG_CHARGED_AA = list('DE')
_AROMATIC = list('HFWY')

_ALL_AA = [_SMALL_HYDROPHOBIC_AA, _LARGE_HYDROPHOBIC_AA, _POLAR_AA,
           _POS_CHARGED_AA, _NEG_CHARGED_AA, _AROMATIC]


def _get_aa_list(aa):
    """Get a list of amino acids similiar to a given amino acid.

    Parameters
    ----------
    aa : str
        Amino acid in one letter code.

    Returns
    -------
    list of str
        List of amino acids similiar to the given one.

    """
    aa_set = set()
    for aa_list in _ALL_AA:
        if aa in aa_list:
            aa_set.update(aa_list)
    return list(aa_set)


class AntibodyDeimmuPreparation(AbstractDeimmuPreparation):

    """Preparation and generation of input files used to de-immunize
    an antibody sequence.

    This class provides two options for setting possible mutations in each
    position: (1) Only backmutating to the corresponding murine amino acid
    is allowed. (2) Amino acids similiar to the murine one are also considered.

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

    def __init__(self, config, alignment, ev_couplings, source_seq,
                 include_similar_amino_acids=False):
        super(AntibodyDeimmuPreparation, self).__init__(config, alignment, ev_couplings)
        self.set_possible_mutations(source_seq, include_similar_amino_acids)
        self.extract_ev_paras(ev_couplings)

    def set_possible_mutations(self, source_seq, include_similar_amino_acids=False):
        """For each wildtype position, generate a set of amino acid
        substitutions (the wildtype amino acid included).

        Method fills the attribute `possible_mutations`.

        Parameters
        ----------
        source_seq : str
            Sequence which serves as source of backmutation residues. Usually,
            this is the murine sequence to be humanized.
        include_similar_amino_acids : bool, optional (default: False)
            True, if amino acids similiar to the `source_seq` residues are
            also considered.

        """
        print 'Possible mutations:'
        mapping = self.wildtype.map_to_seq(source_seq.upper())
        self.possible_mutations = [set(res) for res in self.wildtype.sequence]
        for i, muts in enumerate(self.possible_mutations):
            if i not in self.excluded_pos and i not in self.ignored_pos:
                muts.add(source_seq[mapping[i]])
                if include_similar_amino_acids:
                    similiar_aa = _get_aa_list(source_seq[mapping[i]])
                    muts.update(similiar_aa)
            print muts


def init_params(args):
    """Initialize arguments parsed from the command line."""
    config = ConfigParser.ConfigParser()
    config.read(args.config)
    args.config = config

    args.alignment = AlignIO.read(args.alignment, 'fasta')
    ev_couplings = EVcouplings(args.eij_filename)
    vars(args)['ev_couplings'] = ev_couplings
    del vars(args)['eij_filename']

    source_record = SeqIO.read(args.source_seq, 'fasta')
    args.source_seq = str(source_record.seq).upper()
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
    parser.add_argument('--source_seq', '-s', required=True, help="""Sequence which
                        serves as source of backmutation residues. Usually,
                        this is the murine sequence to be humanized.""")
    parser.add_argument('--out', '-o', required=True, help="""The output file
                        of the data model (in AMPL/GMPL) format.""")
    parser.add_argument('--model', '-m', required=True, nargs=2, help="""The ILP
                        model files in MathProg (first argument: imm, second
                        argument: en).""")
    parser.add_argument('--include_similar_amino_acids', action='store_true',
                        help="""In addition to the murine residues, also consider
                        amino acis similiar to the murine one.""")
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
