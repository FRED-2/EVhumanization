"""Module for representing allels."""

import numpy as np

from EVhumanization.utilities import ALPHABET_PROTEIN_NOGAP
from tools import to_data_format


class AlleleCollection(object):

    """Representation of a set of alleles.

    Parameters
    ----------
    config : `ConfigParser`
        Config file, describing further parameters for de-immunization.

    Attributes
    ----------
    alleles : list of `Allele`
        List of Alleles.

    """

    def __init__(self, config):
        self.read_in_alleles(config.get('sets', 'allele_file'))
        pssm_pos_const = float(config.get('generel', 'pssm_pos_const'))
        epitope_length = int(config.get('parameters', 'epi_len'))
        for allele in reversed(self.alleles):
            is_set = allele.read_in_tepitope_pssm(pssm_pos_const,
                                                  epitope_length)
            if is_set:
                allele.normalize_pssm(pssm_pos_const, epitope_length)
            else:
                self.alleles.remove(allele)
        #print self

    def __str__(self):
        s = '# Allele name\t probability\t pssm_thresh\n'
        for allele in self.alleles:
            s += str(allele) + '\n'
        return s

    def read_in_alleles(self, filename):
        """Reads allele information from file.

        Parameters
        ----------
        filename : str
            Name of the input file.

        """
        self.alleles = []
        with open(filename, 'rU') as allele_file:
            for line in allele_file:
                name, pssm_thresh, probability = line.split(',')
                a = Allele(name, float(probability),
                                    float(pssm_thresh))
                a.read_in_tepitope_pssm(0, 9)
                a.normalize_pssm(0, 9)
                self.alleles.append(a)

    def to_set_A(self):
        """Convert allele names to data format."""
        return ' '.join([allele.name for allele in self.alleles])

    def to_param_pssm_thresh(self):
        """Convert PSSM threshold of the alleles to data format."""
        return '\n' + '\n'.join(
            '\t' + allele.name + '\t' + str(allele.pssm_thresh)
            for allele in self.alleles
        )

    def to_param_p(self):
        """Convert probability of the alleles to data format."""
        return '\n' + '\n'.join(
            '\t' + allele.name + '\t' + str(allele.probability)
            for allele in self.alleles
        )

    def to_param_pssm(self, epitope_length):
        """Convert PSSM of the alleles to data format."""
        return '\n' + ''.join(
            allele.pssm_to_data_format(epitope_length) + '\n'
            for allele in self.alleles
        )


class Allele(object):

    """Representing an allele.

    Parameters
    ----------
    name : str
        name of the allele.
    probability : float
        probability of the allele in the target population.
    pssm_thresh : float
        PSSM threshold.
    pssm : dict (key=tuple of length 2 (str, int), value=float), optional (default: None)
        Position-specific scoring matrix (PSSM),
        as dictionary: key=(aa_i, i), value=PSSM value of amino acid aa_i
        in position i

    Attributes
    ----------
    name : str
        name of the allele.
    probability : float
        probability of the allele in the target population.
    pssm_thresh : float
        PSSM threshold.
    pssm : dict (key=tuple of length 2 (str, int), value=float)
        Position-specific scoring matrix (PSSM),
        as dictionary: key=(aa_i, i), value=PSSM value of amino acid aa_i
        in position i

    """

    def __init__(self, name, probability, pssm_thresh, pssm=None):
        self.name = name
        self.probability = probability
        self.pssm_thresh = pssm_thresh
        self.pssm = pssm

    def __str__(self):
        return 'Allele %s\t%f\t%f' % (self.name, self.probability, self.pssm_thresh)

    def read_in_tepitope_pssm(self, pssm_const, epitope_length):
        """Read TEPITOPE PSSM from file.

        Parameters
        ----------
        pssm_const : float
            PSSM constant.
        epitope_length : int
            Length of an epitope.

        Returns
        -------
        bool
            True, if the reading was successful.

        """
        try:
            name = self.name.replace("HLA-","").replace(":","").replace("*","_")
            pssm_tmp = getattr(__import__("EVhumanization.utilities.tepitopepan_matrices", fromlist=[name]),
                               name)
        except ImportError:
            print 'Allele %s not supported' % self.name
            return False

        self.pssm = {}
        for i in xrange(epitope_length):
            self.pssm[('C', i)] = pssm_const

        for i, val in pssm_tmp.iteritems():
            for aa, v in val.iteritems():
                self.pssm[(aa, i)] = v + pssm_const

        return True

    def normalize_pssm(self, pssm_const, epitope_length):
        """Normalize PSSM by z-Score.

        Parameters
        ----------
        pssm_const : float
            PSSM constant.
        epitope_length : int
            Length of an epitope.

        """
        mean = np.mean(self.pssm.values())
        std = np.std(self.pssm.values())
        self.pssm = {k: (v - mean) / std + pssm_const
                     for k, v in self.pssm.iteritems()}
        self.pssm_thresh = (self.pssm_thresh - epitope_length * mean)\
            / std + epitope_length * pssm_const

    def pssm_to_data_format(self, epitope_length):
        """Convert PSSM of the allele to data format."""
        keyword = '[%s,*,*]:' % self.name
        identifier = ' '.join(str(i + 1) for i in xrange(epitope_length))
        content = '\n' + '\n'.join([
            aa + '\t' + '\t'.join(str(self.pssm[aa, i])
                                  for i in xrange(epitope_length))
                        for aa in list(ALPHABET_PROTEIN_NOGAP)
        ])
        return to_data_format(
            keyword, identifier, content,
            end=False, new_line=False
        )
