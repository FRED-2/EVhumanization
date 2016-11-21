import sys
import subprocess
from abc import ABCMeta, abstractmethod

from alleles import AlleleCollection
from wildtype import Wildtype

sys.path.append('utilities')
from ev_couplings_v4 import ALPHABET_PROTEIN_NOGAP


class AbstractDeimmuPreparation(object):
    """Abstract class for the preparation and generation of input files
    used to de-immunize an amino acid sequence.

    Attributes:
        wildtype -- wildtype to be deimmunized
        num_mutations -- number of allowed mutations
        epitope_length -- length of an epitope
        allele_coll -- collection of alleles
        excluded_pos -- list of positions not allowed to be mutated
        ignored_pos -- list of positions ignored in calculation
        possible_mutations -- list of possible amino acid substitutions for each position
        hi -- his relevant for the calculation,
              as dictionary: key = (position, residue), value = hi value
        eij_indices -- list of eij indices relevant for the calculation
        eij -- eijs relevant for the calculation,
               as dictionary: key = (position 1, position 2, residue 1, residue 2), value = eij value
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
        """
            For each wildtype residue generate a set of possible
            amino acid substitutions.
        """
        raise NotImplementedError

    def read_in_config_paras(self, config):
        """Read in simple parameters from the config file."""
        self.num_mutations = int(config.get('parameters', 'k'))
        self.epitope_length = int(config.get('parameters', 'epi_len'))
        ex_pos = config.get('sets','exclude_pos').split(',')
        ig_pos = config.get('sets','ignore_pos').split(',')
        self.excluded_pos = map(lambda x: int(x.strip())-1, ex_pos) if ex_pos[0] != '' else []
        self.ignored_pos = map(lambda x: int(x.strip())-1, ig_pos) if ig_pos[0] != '' else []

        print 'Number of mutations: param k = %i' % self.num_mutations
        print 'Epitope length: %i' % self.epitope_length
        print 'Excluded positions: %s' % self.excluded_pos
        print 'Ignored positions: %s\n' % self.ignored_pos

    def extract_ev_paras(self, ev_couplings):
        """Extract relevant fields and pair coupling parameters."""
        self.hi = {}
        for i in self.wildtype.epitope_pos:
            for ai in filter(lambda x: x != '-', ev_couplings.alphabet_map.keys()):
                self.hi[(i,ai)] = -1.0 * ev_couplings.h_i[
                    ev_couplings.index_map[self.wildtype.map_to_uniprot(i)],
                    ev_couplings.alphabet_map[ai]
                ]
        self.eij_indices = sorted(set(
            [(i,j) if i < j else (j,i)
             for i in self.wildtype.epitope_pos
             for j in self.wildtype.epitope_pos
             if i != j and ((i < j and i not in self.ignored_pos)
             or (j < i and j not in self.ignored_pos))]
        ))
        self.eij = {}
        for i,j in self.eij_indices:
            for ai in self.possible_mutations[i]:
                for aj in self.possible_mutations[j]:
                    self.eij[(i,j,ai,aj)] = -1.0 * ev_couplings.e_ij[
                        ev_couplings.index_map[self.wildtype.map_to_uniprot(i)],
                        ev_couplings.index_map[self.wildtype.map_to_uniprot(j)],
                        ev_couplings.alphabet_map[ai],
                        ev_couplings.alphabet_map[aj]
                    ]

    def to_data_file(self, out):
        """Generate data file."""
        with open(out, 'w') as f:
            f.write('##################################\n#\n# Sets I\n#\n##################################\n\n')
            f.write('set SIGMA := ' + ' '.join(list(ALPHABET_PROTEIN_NOGAP)) + ';\n')
            f.write('set A := ' + ' '.join([a.name for a in self.allele_coll.alleles]) + ';\n\n')

            f.write('##################################\n#\n# Params I\n#\n##################################\n\n')
            f.write('param N := ' + str(len(self.wildtype.sequence)) + ';\n')
            f.write('param eN := ' + str(self.epitope_length) + ';\n')
            f.write('param k := ' + str(self.num_mutations) + ';\n')
            f.write('param pssm_thresh :=\n' + '\n'.join('\t'+a.name+'\t'+str(a.pssm_thresh) for a in self.allele_coll.alleles) + ';\n\n')
            f.write('param p :=\n' + '\n'.join('\t'+a.name+'\t'+str(a.probability) for a in self.allele_coll.alleles) + ';\n\n')
            f.write('param pssm :=\n')
            for a in self.allele_coll.alleles:
                f.write('[' + a.name + ',*,*]: ' + ' '.join(str(i+1) for i in xrange(self.epitope_length))
                      + ':=\n' + '\n'.join([aa+'\t'+'\t'.join(str(a.pssm[aa,i]) for i in xrange(self.epitope_length)) for aa in list(ALPHABET_PROTEIN_NOGAP)]) + '\n')
            f.write(';\n\n')

            f.write('##################################\n#\n# Sets II\n#\n##################################\n\n')
            f.write('set Eij := ' + '\t'.join(str(i+1) + ' ' + str(j+1) for i,j in self.eij_indices) + ';\n')
            f.write('set E := '   + ' '.join(str(i+1) for i in self.wildtype.epitope_pos) + ';\n\n')
            for i in xrange(len(self.wildtype.sequence)):
                f.write('set WT[' + str(i+1) + '] := ' + self.wildtype.sequence[i] + ';\n')
            f.write('\n'.join('set M[' + str(i+1) + '] := ' + ' '.join(a for a in self.possible_mutations[i]) + ';'\
                for i in xrange(len(self.wildtype.sequence))) + '\n\n')

            f.write('##################################\n#\n# Params II\n#\n##################################\n\n')
            f.write('param h: '+' '.join(a for a in list(ALPHABET_PROTEIN_NOGAP))+' := \n')
            f.write('\n'.join(str(i+1) + ' ' + ' '.join(str(self.hi[i,a])\
                for a in list(ALPHABET_PROTEIN_NOGAP)) for i in self.wildtype.epitope_pos)+';\n\n')
            f.write('param eij :=\n')
            for i,j in self.eij_indices:
                f.write('[' + str(i+1) + ',' + str(j+1) + ',*,*]: ' + '\t'.join(self.possible_mutations[j]) + ' :=\n')
                for ai in self.possible_mutations[i]:
                    f.write(ai + '\t' + '\t'.join(str(self.eij[i,j,ai,aj]) for aj in self.possible_mutations[j]) + '\n')
            f.write(';\n\n')
            f.write('end;\n\n')

    def generate_lp_files(self, out, model):
        """Generate input files for the solver."""
        lp_model = ''.join(out.split('.')[:-1])
        print 'glpsol -m %s -d %s --check --wlp %s'%(model,'.'.join(out.split('.')[:-1])+'_imm.data',lp_model)
        subprocess.call('glpsol -m %s -d %s --check --wcpxlp %s'%(model[0],out,lp_model+'_imm.lp'), shell=True) # Benni: wlp
        subprocess.call('glpsol -m %s -d %s --check --wcpxlp %s'%(model[1],out,lp_model+'_en.lp'), shell=True)
