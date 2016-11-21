import os
import numpy as np
from abc import ABCMeta, abstractmethod


class AlleleCollection(object):

    def __init__(self, config):
        self.read_in_alleles(config.get('sets', 'allele_file'))
        pssm_dir = config.get('parameters', 'pssm_dir')
        pssm_pos_const = float(config.get('generel', 'pssm_pos_const'))
        epitope_length = int(config.get('parameters', 'epi_len'))
        for allele in reversed(self.alleles):
            is_set = allele.read_in_tepitope_pssm(pssm_dir, pssm_pos_const, epitope_length)
            if is_set:
                allele.normalize_pssm(pssm_pos_const, epitope_length)
            else:
                self.alleles.remove(allele)
        self.print_info()

    def read_in_alleles(self, filename):
        """Reads allele information from file."""
        self.alleles = []
        with open(filename, 'rU') as allele_file:
            for line in allele_file:
                name, pssm_thresh, probability = line.split(',')
                self.alleles.append(Allele(name, float(probability), float(pssm_thresh)))

    def print_info(self):
        print '# Allele name\t probability\t pssm_thresh'
        for allele in self.alleles:
            print allele
        print


class Allele(object):
    """Representing an allele.

    attributes:
        name -- name of the allele
        probability -- probability of the allele in the target population
        pssm_thresh -- TODO
        pssm -- position-specific scoring matrix
    """

    def __init__(self, name, probability, pssm_thresh, pssm=None):
        self.name = name
        self.probability = probability
        self.pssm_thresh = pssm_thresh
        self.pssm = pssm

    def __str__(self):
        return 'Allele %s\t%f\t%f' % (self.name, self.probability, self.pssm_thresh)

    def read_in_tepitope_pssm(self, pssms_dir, pssm_const, epitope_length):
        """
            Reads TEPITOPE pssm from file and returns wether the reading
            was successful or not.
        """
        tepitope_aa_list = list('ADEFGHIKLMNPQRSTVWY')
        if not os.path.exists(pssms_dir + self.name):
            print 'Allele %s not supported' % self.name
            return False
        self.pssm = {}
        for i in xrange(epitope_length):
            self.pssm[('C',i)] = pssm_const
        with open(pssms_dir + self.name, 'rU') as pssm_file:
            for aa_index, line in enumerate(pssm_file):
                pssm_values = line.split()
                for i in xrange(epitope_length):
                    self.pssm[(tepitope_aa_list[aa_index],i)] = float(pssm_values[i]) + pssm_const
        return True

    def normalize_pssm(self, pssm_const, epitope_length):
        """Normalizes pssm by Z-Score."""
        mean = np.mean(self.pssm.values())
        std = np.std(self.pssm.values())
        self.pssm = {k:(v-mean)/std+pssm_const for k, v in self.pssm.iteritems()}
        self.pssm_thresh = (self.pssm_thresh-epitope_length*mean)/std + epitope_length*pssm_const
