"""Base class for the preparation and generation of input files used to
de-immunize an amino acid sequence.

Any subclass needs to implement the abstract method `set_possible_mutations`,
filling the attribute `possible_mutations`.

"""

import sys
import subprocess
from abc import ABCMeta, abstractmethod

from alleles import AlleleCollection
from wildtype import Wildtype
from utilities.ev_couplings_v4 import ALPHABET_PROTEIN_NOGAP
from tools import to_data_format


class AbstractDeimmuPreparation(object):

    """Abstract class for the preparation and generation of input files
    used to de-immunize an amino acid sequence.

    Parameters
    ----------
    config : `ConfigParser`
        Config file, describing further parameters of the de-immunization process.
    alignment : `MultipleSeqAlignment`
        Alignment which was used for plmc calculation in focus mode
        (with upper and lower letters).
    ev_couplings : `EVcouplings`
        e_ij binaries read in, generated from `alignment`.

    Attributes
    ----------
    wildtype : `Wildtype`
        The wildtype to be de-immunized.
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
        as dictionary: key=(i, aa_i), value=h_i of amino acid aa_i in position i.
    eij_indices : list of tuple of length 2 (int, int)
        List of paired e_ij indices whose associated values are relevant
        to the calculation.
    eij : dict (key=tuple of length 4 (int, int, str, str), value=float)
        e_ij pair couplings relevant to the calculation,
        as dictionary: key=(i, j, aa_i, aa_j), value=e_ij of amino acids
        aa_i and aa_j in positions i and j, respectively.

    Methods
    -------
    set_possible_mutations()
        For each wildtype position, generate a set of amino acid
        substitutions (the wildtype amino acid included).

    """

    __metaclass__ = ABCMeta

    def __init__(self, config, alignment, ev_couplings):
        self.read_in_config_paras(config)
        self.allele_coll = AlleleCollection(config)
        self.wildtype = Wildtype.create(alignment, ev_couplings)
        self.wildtype.predict_epitope_positions(self.allele_coll, self.epitope_length)
        self.wildtype.print_epitope_prediction(self.excluded_pos, self.ignored_pos)

    @abstractmethod
    def set_possible_mutations(self):
        """For each wildtype position, generate a set of amino acid
        substitutions (the wildtype amino acid included).

        Fill the attribute `possible_mutations`. Needs to be implemented
        by any subclass.

        Raises
        ------
        NotImplementedError
            If the method is not implemented.

        """
        raise NotImplementedError

    def read_in_config_paras(self, config):
        """Read in simple parameters from the config file.

        Parameters
        ----------
        config : `ConfigParser`
            Config file, describing further parameters for de-immunization.

        """
        self.num_mutations = int(config.get('parameters', 'k'))
        self.epitope_length = int(config.get('parameters', 'epi_len'))
        ex_pos = config.get('sets', 'exclude_pos').split(',')
        ig_pos = config.get('sets', 'ignore_pos').split(',')
        self.excluded_pos = map(lambda x: int(x.strip()) - 1, ex_pos)\
            if ex_pos[0] != '' else []
        self.ignored_pos = map(lambda x: int(x.strip()) - 1, ig_pos)\
            if ig_pos[0] != '' else []

        print 'Number of mutations: param k = %i' % self.num_mutations
        print 'Epitope length: %i' % self.epitope_length
        print 'Excluded positions: %s' % self.excluded_pos
        print 'Ignored positions: %s\n' % self.ignored_pos

    def extract_ev_paras(self, ev_couplings):
        """Extract relevant h_i fields and e_ij pair couplings.

        Parameters
        ----------
        ev_couplings : `EVcouplings`
            e_ij binaries read in, generated from `alignment`.

        """
        self.hi = {}
        for i in self.wildtype.epitope_pos:
            for ai in filter(lambda x: x != '-', ev_couplings.alphabet_map.keys()):
                self.hi[(i, ai)] = -1.0 * ev_couplings.h_i[
                    ev_couplings.index_map[self.wildtype.map_to_uniprot(i)],
                    ev_couplings.alphabet_map[ai]
                ]
        self.eij_indices = sorted(set(
            [(i, j) if i < j else (j, i)
             for i in self.wildtype.epitope_pos
             for j in self.wildtype.epitope_pos
             if i != j and ((i < j and i not in self.ignored_pos) or
                            (j < i and j not in self.ignored_pos))]
        ))
        self.eij = {}
        for i, j in self.eij_indices:
            for ai in self.possible_mutations[i]:
                for aj in self.possible_mutations[j]:
                    self.eij[(i, j, ai, aj)] = -1.0 * ev_couplings.e_ij[
                        ev_couplings.index_map[self.wildtype.map_to_uniprot(i)],
                        ev_couplings.index_map[self.wildtype.map_to_uniprot(j)],
                        ev_couplings.alphabet_map[ai],
                        ev_couplings.alphabet_map[aj]
                    ]

    def to_data_file(self, out):
        """Generate data file.

        Parameters
        ----------
        out : str
            Name of the output file.

        """
        with open(out, 'w') as f:
            f.write('##################################\n#\n# Sets I\n#\n##################################\n\n')
            f.write(to_data_format('set', 'SIGMA', ' '.join(list(ALPHABET_PROTEIN_NOGAP))))
            f.write(to_data_format('set', 'A', self.allele_coll.to_set_A()) + '\n')

            f.write('##################################\n#\n# Params I\n#\n##################################\n\n')
            f.write(to_data_format('param', 'N', self.wildtype.to_param_N()))
            f.write(to_data_format('param', 'eN', str(self.epitope_length)))
            f.write(to_data_format('param', 'k', str(self.num_mutations)))
            f.write(to_data_format('param', 'pssm_thresh', self.allele_coll.to_param_pssm_thresh()) + '\n')
            f.write(to_data_format('param', 'p', self.allele_coll.to_param_p()) + '\n')
            f.write(to_data_format('param', 'pssm', self.allele_coll.to_param_pssm(self.epitope_length)) + '\n')

            f.write('##################################\n#\n# Sets II\n#\n##################################\n\n')
            f.write(to_data_format('set', 'Eij', self.format_eij_indices()))
            f.write(to_data_format('set', 'E', self.wildtype.to_set_E()) + '\n')
            for i in xrange(len(self.wildtype.sequence)):
                f.write(to_data_format('set', 'WT[%i]' % (i + 1), self.wildtype.sequence[i]))
            for i in xrange(len(self.wildtype.sequence)):
                f.write(to_data_format('set', 'M[%i]' % (i + 1), ' '.join(aa for aa in self.possible_mutations[i])))
            f.write('\n')

            f.write('##################################\n#\n# Params II\n#\n##################################\n\n')
            f.write(to_data_format('param', 'h' + ': ' + ' '.join(aa for aa in list(ALPHABET_PROTEIN_NOGAP)), self.format_hi()) + '\n')
            f.write(to_data_format('param', 'eij', self.format_eij()))
            f.write('\nend;\n\n')

    def format_eij_indices(self):
        """Convert e_ij indices to data format."""
        return '\t'.join(str(i + 1) + ' ' + str(j + 1)
                         for i, j in self.eij_indices)

    def format_hi(self):
        """Convert h_i fields to data format."""
        return '\n' + '\n'.join(str(i + 1) + ' ' + ' '.join(
            str(self.hi[i, aa])
            for aa in list(ALPHABET_PROTEIN_NOGAP))
            for i in self.wildtype.epitope_pos
        )

    def format_eij(self):
        """Convert e_ij pair couplings to data format."""
        s = '\n'
        for i, j in self.eij_indices:
            keyword = '[%i,%i,*,*]:' % (i + 1, j + 1)
            identifier = '\t'.join(self.possible_mutations[j])
            content = '\n'
            for ai in self.possible_mutations[i]:
                content += ai + '\t' + '\t'.join(
                    str(self.eij[i, j, ai, aj])
                    for aj in self.possible_mutations[j]
                ) + '\n'
            s += to_data_format(keyword, identifier, content, end=False, new_line=False)
        return s

    def generate_lp_files(self, out, model):
        """Generate input files for the solver.

        Parameters
        ----------
        out : str
            Name of the ouput file.
        model : list of str
            Model files.

        """
        lp_model = ''.join(out.split('.')[:-1])
        print ('glpsol -m %s -d %s --check --wlp %s'
               % (model, '.'.join(out.split('.')[:-1]) + '_imm.data', lp_model))
        subprocess.call('glpsol -m %s -d %s --check --wcpxlp %s'
                        % (model[0], out, lp_model + '_imm.lp'), shell=True)  # Benni: wlp
        subprocess.call('glpsol -m %s -d %s --check --wcpxlp %s'
                        % (model[1], out, lp_model + '_en.lp'), shell=True)
