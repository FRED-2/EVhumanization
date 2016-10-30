#!/usr/local/bin/python2.7
# encoding: utf-8

"""
Deimmunize -- shortdesc

Deimmunize is a description

It defines classes_and_methods

@author:     user_name
@copyright:  2013 organization_name. All rights reserved.
@license:    license
@contact:    user_email
@deffield    updated: Updated
"""

import os
import argparse
import ConfigParser
import numpy as np
import subprocess
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

from ev_couplings import EVcouplings, AA_MAP, AA_LIST


def predictstart(x):
    return predict(*x)

def predict(func, args):
    return func(*args)

def tepitope_pred(pos, ep, pssm, allele, thr):
    is_binder = sum(pssm[allele, ep[i], i] for i in xrange(len(ep))) >= thr
    return (allele,pos) if is_binder else (allele,-1)


class Allele(object):
    """Representing an allele.

    attributes:
        name -- name of the allele
        probability -- probability of the allele in the target population
        pssm_thresh -- TODO
        pssm -- position-specific scoring matrix
    """

    _TEPITOPE_AA_LIST = ['A','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    def __init__(self, name, probability, pssm_thresh, pssm=None):
        self.name = name
        self.probability = probability
        self.pssm_thresh = pssm_thresh
        self.pssm = pssm

    def __str__(self):
        return 'Allele %s\t%f\t%f' % (self.name, self.probability, self.pssm_thresh)

    # TODO generalize (specified for TEPITOPE right now)
    def set_pssm(self, pssms_dir, pssm_const, epitope_length):
        """Reads pssm from file and returns wether the reading was successful or not."""
        if not os.path.exists(pssms_dir + self.name):
            print 'Allele %s not supported' % self.name
            return False
        pssm_const = 0
        self.pssm = {}
        for i in xrange(epitope_length):
            self.pssm[('C',i)] = pssm_const
        with open(pssms_dir + self.name, 'rU') as pssm_file:
            for aa_index, line in enumerate(pssm_file):
                pssm_values = line.split()
                for i in xrange(epitope_length):
                    self.pssm[(Allele._TEPITOPE_AA_LIST[aa_index],i)] = float(pssm_values[i]) + pssm_const
        return True

    def normalize_pssm(self, pssm_const, epitope_length):
        """Normalizes pssm by Z-Score."""
        mean = np.mean(self.pssm.values())
        std = np.std(self.pssm.values())
        self.pssm = {k:(v-mean)/std+pssm_const for k, v in self.pssm.iteritems()}
        self.pssm_thresh = (self.pssm_thresh-epitope_length*mean)/std + epitope_length*pssm_const


class Wildtype(object):
    """Representing the wildtype to be deimmunized.

    attributes:
        name -- name of the wildtype
        sequence -- sequence of the wildtype
        index_to_uniprot_offset -- TODO
        uniprot_offset -- TODO
        internal_offset -- TODO
        num_epitopes -- number of predicted epitopes
        epitope_pos -- set of positions belonging to an epitope
        epitope_pos_alleles -- positions belonging to an epitope with specified allele,
                               as dictionary: key = name of the allele, value = set of positions
    """

    def __init__(self, name, sequence, index_to_uniprot_offset=None, uniprot_offset=None,
                 internal_offset=None, num_epitopes=None, epitope_pos=None, epitope_pos_alleles=None):
        self.name = name
        self.sequence = sequence
        self.index_to_uniprot_offset = index_to_uniprot_offset
        self.uniprot_offset = uniprot_offset
        self.internal_offset = internal_offset
        self.num_epitopes = num_epitopes
        self.epitope_pos = epitope_pos
        self.epitope_pos_alleles = epitope_pos_alleles

    def __str__(self):
        s = '# Wildtype to be deimmunized\n'
        s += '>' + self.name + '\n' + self.sequence
        s += '\nLength: %i' % len(self.sequence)
        s += '\nOffsets: (uniprot=%i, internal=%i)' % (self.uniprot_offset, self.internal_offset)
        return s

    def calc_offsets(self, ev_couplings):
        """Calculates uniprot and internal offset of the wildtype sequence."""
        # uniprot offset
        self.index_to_uniprot_offset = ev_couplings.index_to_uniprot_offset
        self.uniprot_offset = self.index_to_uniprot_offset[0]
        # internal offset
        self.internal_offset = 0
        for i in xrange(len(self.sequence)):
            if self.sequence[i].isupper():
                self.internal_offset = i
                break
        self.sequence = self.sequence.upper()

    def predict_epitope_positions(self, alleles, epitope_length):
        """Predicts epitope positions within the wildtype sequence."""
        pssms = {(allele.name, k[0], k[1]):v for allele in alleles for k, v in allele.pssm.iteritems()}

        # get starting positions of epitopes (with specified alleles)
        pool = Pool(cpu_count())
        epitope_start_pos_alleles = pool.map(predictstart, [(tepitope_pred,(i,self.sequence[i:i+9],pssms,allele.name,allele.pssm_thresh))\
            for i in xrange(len(self.sequence) - (epitope_length-1))\
            for allele in alleles])
        self.num_epitopes = len(filter(lambda x: x[1] != -1, epitope_start_pos_alleles))

        self.epitope_pos_alleles = defaultdict(set)
        for allele_name, pos in epitope_start_pos_alleles:
            if pos > 0:
                for j in range(pos,pos+9):
                    self.epitope_pos_alleles[allele_name].add(j)

        # epitope_pos_ext is the extended epitope set with buffer of 8 aa at the beginning of each epitope
        epitope_start_pos = set([pos for _,pos in epitope_start_pos_alleles])
        self.epitope_pos, epitope_pos_ext = [], []
        for i in epitope_start_pos:
            if i >= 0:
                self.epitope_pos.extend(range(i,i+epitope_length))
                if i >= (epitope_length-2):
                    end = i+epitope_length if i+epitope_length < (len(self.sequence)-8) else len(self.sequence)-8
                    epitope_pos_ext.extend(range(i-(epitope_length-2),end))
                else:
                    epitope_pos_ext.extend(range(i+epitope_length))
        self.epitope_pos.sort(), epitope_pos_ext.sort()
        self.epitope_pos, epitope_pos_ext = set(self.epitope_pos), set(epitope_pos_ext)

        # map EVFold uniprot index to internal index
        epitope_pos_mapped = set(self.map_to_seq_idx(pos) for pos in self.index_to_uniprot_offset)
        if len(epitope_pos_mapped & self.epitope_pos) == 0:
            print >> stderr, "No epitopes in area considered during Eij calculation! Apply normal epitope reduction."
            exit(1)
        self.epitope_pos = epitope_pos_mapped

    def print_epitope_prediction(self, excluded_pos, ignored_pos):
        """Prints epitope prediction of the wildtype sequence for each allele."""
        s = 'Epitope Prediction of the wildtype\n'
        s += 'Number of predictet epitopes: %i\n' % self.num_epitopes
        s += str(self.epitope_pos) + '\n'
        for name, pos in self.epitope_pos_alleles.iteritems():
            s += name + ': ' + str(pos) + '\n'
        for allele_name, epitope_pos in self.epitope_pos_alleles.iteritems():
            s += allele_name + '\t:'
            for i in range(len(self.sequence)):
                if i in epitope_pos:
                    s += 'x' if i in excluded_pos or i in ignored_pos else '*'
                elif i in excluded_pos:
    		        s += 'x'
                else:
                    s += '-'
    	    s += '\n'
        print s

    def map_to_seq_idx(self, uniprot_idx):
        return self.internal_offset + (uniprot_idx - self.uniprot_offset)

    def map_to_uniprot(self, pos):
        return self.uniprot_offset + (pos - self.internal_offset)


class Deimmunization(object):
    """Preparation and generation of input files for the de-immunization of an amino acid sequence.

    attributes:
        wildtype -- wildtype to be deimmunized
        num_mutations -- number of allowed mutations
        epitope_length -- length of an epitope
        alleles -- list of alleles
        excluded_pos -- list of positions not allowed to be mutated
        ignored_pos -- list of positions ignored in calculation
        possible_mutations -- list of possible amino acid substitutions for each position
        hi -- his relevant for the calculation,
              as dictionary: key = (position, residue), value = hi value
        eij_indices -- list of eij indices relevant for the calculation
        eij -- eijs relevant for the calculation,
               as dictionary: key = (position 1, position 2, residue 1, residue 2), value = eij value
    """

    def __init__(self, args):

        config = ConfigParser.ConfigParser()
        config.read(args.config)

        alignment = AlignIO.read(args.alignment, 'fasta')

        ev_couplings = EVcouplings(args.eij)
        ev_couplings.to_zero_sum_gauge()

        # get wildtype to be deimmunized as first record in the given alignment
        self.wildtype = Wildtype(alignment[0].id, str(alignment[0].seq))
        self.wildtype.calc_offsets(ev_couplings)

        self.num_mutations = int(config.get('parameters', 'k'))
        self.epitope_length = int(config.get('parameters', 'epi_len'))
        self.alleles = Deimmunization.read_in_alleles(config.get('sets', 'allele_file'))

        print str(self.wildtype) + '\n'
        print 'Number of mutations: param k = %i\nEpitope length: %i\n' % (self.num_mutations, self.epitope_length)
        print '# Allele name\t probability\t pssm_thresh'
        for allele in self.alleles:
            print allele
        print

        # get pssm for each allele
        for allele in reversed(self.alleles):
            is_set = allele.set_pssm(config.get('parameters', 'pssm_dir'), float(config.get('generel', 'pssm_pos_const')), self.epitope_length)
            if is_set:
                allele.normalize_pssm(float(config.get('generel', 'pssm_pos_const')), self.epitope_length)
            else:
                self.alleles.remove(allele)

        # get excluded and ignored positions
        ex_pos = config.get('sets','exclude_pos').split(',')
        ig_pos = config.get('sets','ignore_pos').split(',')
        self.excluded_pos = map(lambda x: int(x.strip())-1, ex_pos) if ex_pos[0] != '' else []
        self.ignored_pos = map(lambda x: int(x.strip())-1, ig_pos) if ig_pos[0] != '' else []
        print 'Excluded positions: %s\nIgnored positions: %s\n' % (self.excluded_pos, self.ignored_pos)

        # predict epitope positions within the wildtype sequence
        self.wildtype.predict_epitope_positions(self.alleles, self.epitope_length)
        self.wildtype.print_epitope_prediction(self.excluded_pos, self.ignored_pos)

        # create list of amino acids considered for mutations in each position
        print 'Possible mutations:'
        self.set_possible_mutations(alignment, float(config.get('generel', 'frequency_thresh')))

        # extract relevant his
        self.hi = {}
        for i in self.wildtype.epitope_pos:
            for ai in filter(lambda x: x != '-', AA_MAP.keys()):
                self.hi[(i,ai)] = -1.0*ev_couplings.fields[ev_couplings.map_uniprot_index[self.wildtype.map_to_uniprot(i)], AA_MAP[ai]]

        # extract relevant eij indices
        self.eij_indices = sorted(set([(i,j) if i < j else (j,i)\
                               for i in self.wildtype.epitope_pos for j in self.wildtype.epitope_pos\
                                   if i != j and ((i < j and i not in self.ignored_pos) or (j < i and j not in self.ignored_pos))]))

        # extract relevant eijs
        self.eij = {}
        for i,j in self.eij_indices:
            for ai in self.possible_mutations[i]:
                for aj in self.possible_mutations[j]:
                    self.eij[(i,j,ai,aj)] = -1.0*ev_couplings.pair_couplings[ev_couplings.map_uniprot_index[self.wildtype.map_to_uniprot(i)],\
                                                                             ev_couplings.map_uniprot_index[self.wildtype.map_to_uniprot(j)],\
                                                                             AA_MAP[ai], AA_MAP[aj]]


    @staticmethod
    def read_in_alleles(filename):
        """Reads allele information from file."""
        alleles = []
        with open(filename, 'rU') as allele_file:
            for line in allele_file:
                name, pssm_thresh, probability = line.split(',')
                alleles.append(Allele(name, float(probability), float(pssm_thresh)))
        return alleles

    def set_possible_mutations(self, alignment, freq_thresh):
        self.possible_mutations = [None]*len(self.wildtype.sequence)
        for i in set(range(len(self.wildtype.sequence))) - self.wildtype.epitope_pos:
            self.possible_mutations[i] = set(self.wildtype.sequence[i])

        aln_summ = SummaryInfo(alignment)
        all_letters = aln_summ._get_all_letters()
        for char in ['.','-']:
            all_letters = all_letters.replace(char, '')
        for j in self.wildtype.epitope_pos:
            freq = aln_summ._get_letter_freqs(j, aln_summ.alignment._records, all_letters, ['.','-'])
            tmp = set()
            if j not in self.excluded_pos and j not in self.ignored_pos:
                tmp = set([aa.upper() for aa,fr in freq.items() if fr > freq_thresh and aa not in ['-','.','X','x']])
            if self.wildtype.sequence[j] not in tmp:
                tmp.add(self.wildtype.sequence[j])
            self.possible_mutations[j] = tmp
            print self.possible_mutations[j]

    def to_data_file(self, out, model):
        """Generates data file and input files for the solver."""
        with open(out, 'w') as f:
            f.write('##################################\n#\n# Sets I\n#\n##################################\n\n')
            f.write('set SIGMA := ' + ' '.join(AA_LIST[1:]) + ';\n')
            f.write('set A := ' + ' '.join([a.name for a in self.alleles]) + ';\n\n')

            f.write('##################################\n#\n# Params I\n#\n##################################\n\n')
            f.write('param N := ' + str(len(self.wildtype.sequence)) + ';\n')
            f.write('param eN := ' + str(self.epitope_length) + ';\n')
            f.write('param k := ' + str(self.num_mutations) + ';\n')
            f.write('param pssm_thresh :=\n' + '\n'.join('\t'+a.name+'\t'+str(a.pssm_thresh) for a in self.alleles) + ';\n\n')
            f.write('param p :=\n' + '\n'.join('\t'+a.name+'\t'+str(a.probability) for a in self.alleles) + ';\n\n')
            f.write('param pssm :=\n')
            for a in self.alleles:
                f.write('[' + a.name + ',*,*]: ' + ' '.join(str(i+1) for i in xrange(self.epitope_length))
                      + ':=\n' + '\n'.join([aa+'\t'+'\t'.join(str(a.pssm[aa,i]) for i in xrange(self.epitope_length)) for aa in AA_LIST[1:]]) + '\n')
            f.write(';\n\n')

            f.write('##################################\n#\n# Sets II\n#\n##################################\n\n')
            f.write('set Eij := ' + '\t'.join(str(i+1) + ' ' + str(j+1) for i,j in self.eij_indices) + ';\n')
            f.write('set E := '   + ' '.join(str(i+1) for i in self.wildtype.epitope_pos) + ';\n\n')
            for i in xrange(len(self.wildtype.sequence)):
                f.write('set WT[' + str(i+1) + '] := ' + self.wildtype.sequence[i] + ';\n')
            f.write('\n'.join('set M[' + str(i+1) + '] := ' + ' '.join(a for a in self.possible_mutations[i]) + ';'\
                for i in xrange(len(self.wildtype.sequence))) + '\n\n')

            f.write('##################################\n#\n# Params II\n#\n##################################\n\n')
            f.write('param h: '+' '.join(a for a in AA_LIST[1:])+' := \n')
            f.write('\n'.join(str(i+1) + ' ' + ' '.join(str(self.hi[i,a])\
                for a in AA_LIST[1:]) for i in self.wildtype.epitope_pos)+';\n\n')
            f.write('param eij :=\n')
            for i,j in self.eij_indices:
                f.write('[' + str(i+1) + ',' + str(j+1) + ',*,*]: ' + '\t'.join(self.possible_mutations[j]) + ' :=\n')
                for ai in self.possible_mutations[i]:
                    f.write(ai + '\t' + '\t'.join(str(self.eij[i,j,ai,aj]) for aj in self.possible_mutations[j]) + '\n')
            f.write(';\n\n')
            f.write('end;\n\n')

        # generate input files for the solver
        lp_model = ''.join(out.split('.')[:-1])
        print 'glpsol -m %s -d %s --check --wlp %s'%(model,'.'.join(out.split('.')[:-1])+'_imm.data',lp_model)
        subprocess.call('glpsol -m %s -d %s --check --wcpxlp %s'%(model[0],out,lp_model+'_imm.lp'), shell=True) # Benni: wlp
        subprocess.call('glpsol -m %s -d %s --check --wcpxlp %s'%(model[1],out,lp_model+'_en.lp'), shell=True)


def command_line():
    parser = argparse.ArgumentParser(description='Deimmunization')
    parser.add_argument('--config','-c',
                      required=True,
                      help='Config File')
    parser.add_argument('--eij','-e',
                      required=True,
                      help='Eij File')
    parser.add_argument('--alignment','-a',
                      required=True,
                      help='The complete alignment generated by EVFold with upper and lower letters. '
                          +'The sequence to be deimmunized must be the first in the alignment.')
    parser.add_argument('--out','-o',
                      required=True,
                      help='The output file of the data model (in AMPL/GMPL) format')
    parser.add_argument('--model','-m',
                      required=True,
                      nargs=2,
                      help='The ILP model files in MathProg (first argument: imm, second argument: en)')
    args = parser.parse_args()

    deimmu = Deimmunization(args)
    deimmu.to_data_file(args.out, args.model)

if __name__ == '__main__':
    command_line()
