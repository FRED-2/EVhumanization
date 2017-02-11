#!/usr/local/bin/python2.7
# encoding: utf-8
"""Preparation and generation of input files used to de-immunize an amino
acid sequence.

"""

import sys
from argparse import ArgumentParser
import ConfigParser
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

from utilities.ev_couplings_v4 import EVcouplings
from abstract_deimmu_prep import AbstractDeimmuPreparation


class DeimmuPreparation(AbstractDeimmuPreparation):

    """Preparation and generation of input files used to de-immunize an amino
    acid sequence.

    Parameters
    ----------
    config : `ConfigParser`
        Config file, describing further parameters for de-immunization.
    alignment : `MultipleSeqAlignment`
        Alignment, used for plmc calculation in focus mode
        (with upper and lower letters).
    ev_couplings : `EVcouplings`
        e_ij binaries read in, generated from `alignment`.

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

    def __init__(self, config, alignment, ev_couplings):
        super(DeimmuPreparation, self).__init__(config, alignment, ev_couplings)
        self.set_possible_mutations(alignment, float(config.get('generel', 'frequency_thresh')))
        self.extract_ev_paras(ev_couplings)

    def set_possible_mutations(self, alignment, freq_thresh):
        """For each wildtype position, generate a set of amino acid
        substitutions (the wildtype amino acid included).

        Method fills the attribute `possible_mutations`.

        Parameters
        ----------
        alignment : `MultipleSeqAlignment`
            Alignment, used for plmc calculation in focus mode
            (with upper and lower letters).
        freq_thresh : float
            If the the frequency of amino acid a_i in position i is greater
            than the given frequency threshold, then the amino acid is
            considered as possible mutation in position i.

        """
        print 'Possible mutations:'

        self.possible_mutations = [None] * len(self.wildtype.sequence)
        for i in set(range(len(self.wildtype.sequence))) - self.wildtype.epitope_pos:
            self.possible_mutations[i] = set(self.wildtype.sequence[i])

        aln_summ = SummaryInfo(alignment)
        all_letters = aln_summ._get_all_letters()
        for char in ['.', '-']:
            all_letters = all_letters.replace(char, '')
        for j in self.wildtype.epitope_pos:
            freq = aln_summ._get_letter_freqs(j, aln_summ.alignment._records,
                                              all_letters, ['.', '-'])
            tmp = set()
            if j not in self.excluded_pos and j not in self.ignored_pos:
                tmp = set([aa.upper() for aa, fr in freq.items()
                           if fr > freq_thresh and aa not in ['-', '.', 'X', 'x']])
            if self.wildtype.sequence[j] not in tmp:
                tmp.add(self.wildtype.sequence[j])
            self.possible_mutations[j] = tmp
            print self.possible_mutations[j]


def init_params(args):
    """Initialize arguments parsed from the command line."""
    config = ConfigParser.ConfigParser()
    config.read(args.config)
    args.config = config

    args.alignment = AlignIO.read(args.alignment, 'fasta')

    ev_couplings = EVcouplings(args.eij_filename)
    vars(args)['ev_couplings'] = ev_couplings
    del vars(args)['eij_filename']
    return args


def command_line():
    """Define parser to read arguments from the command line."""
    parser = ArgumentParser(description="""Preparation and generation
                            of input files used to de-immunize an amino acid
                            sequence.""")
    parser.add_argument('--config', '-c', required=True, help="""Config file,
                        describing further parameters for de-immunization.""")
    parser.add_argument('--eij_filename', '-e', required=True, help='eij binary file.')
    parser.add_argument('--alignment', '-a', required=True, help="""Alignment,
                        used for plmc calculation in focus mode (with upper and
                        lower letters). The sequence to be de-immunized must be
                        the first record in the alignment.""")
    parser.add_argument('--out', '-o', required=True, help="""The output file
                        of the data model (in AMPL/GMPL) format.""")
    parser.add_argument('--model', '-m', required=True, nargs=2, help="""The ILP
                        model files in MathProg (first argument: imm, second
                        argument: en).""")
    return init_params(parser.parse_args())


def main(args):
    """Call amino acid de-immunization preparation and write data and lp files."""
    sub_args = {k: v for k, v in vars(args).iteritems() if k not in ['model', 'out']}
    deimmu = DeimmuPreparation(**sub_args)
    deimmu.to_data_file(args.out)
    deimmu.generate_lp_files(args.out, args.model)


if __name__ == '__main__':
    args = command_line()
    main(args)
