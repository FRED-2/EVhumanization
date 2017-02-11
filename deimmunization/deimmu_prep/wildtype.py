"""Module for representing the wildtype to be de-immunized.

This module also provides the functionality to predict epitopes
of the wildtype sequence.

"""

import sys
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62


def _predictstart(x):
    return _predict(*x)


def _predict(func, args):
    return func(*args)


def _tepitope_pred(position, epitope, pssms, allele_name, threshold):
    """TEPITOPE prediction of binding of an epitope to an allele.

    Parameters
    ----------
    position : int
        Starting position of epitope in sequence.
    epitope: str
        The epitope.
    pssms : dict (key=tuple of length 3 (str, str, int), value=float)
        PSSMs as dictionary: key=(allele name, aa_i, i), value=PSSM value
        regarding the specified allele of amino acid aa_i in position i.
    allele_name : str
        Name of the allele.
    threshold : float
        PSSM score threshold to decide if an epitope is a binder or non-binder.

    Returns
    -------
    tuple of length 3 (str, int, float)
        Specification of the PSSM score of an epitope starting at a certain
        position, as tuple: (allele name, starting position of the epitope
        in the sequence, PSSM score)

    """
    pssm_score = sum(pssms[allele_name, epitope[i], i] for i in xrange(len(epitope)))
    return (allele_name, position, pssm_score - threshold)\
        if pssm_score >= threshold else (allele_name, -1, pssm_score)


class Wildtype(object):

    """Representing the wildtype to be de-immunized.

    Parameters
    ----------
    name : str
        Name of the wildtype.
    sequence : str
        Sequence of the wildtype.

    Attributes
    ----------
    name : str
        Name of the wildtype.
    sequence : str
        Sequence of the wildtype.
    index_list : list of int, optional (default: None)
        EVcouplings indices of wildtype sequence.
    uniprot_offset : int, optional (default: None)
        Uniprot offset.
    internal_offset : int, optional (default: None)
        Internal offset.
    num_epitopes : int, optional (default: None)
        Number of predicted epitopes.
    epitope_start_pos_alleles : list of tuple (str, int, float)
        For each allele, lists PSSM score of an epitope defined by its starting
        position in the sequence. The tuples of the list contain the following
        information: (allele name, starting position of an epitope, PSSM score).
    epitope_pos : set of int
        All positions that are part of an epitope.
    epitope_pos_alleles : dict (key=str, value=set of int)
        For each allele, specifies the positions that are part of an epitope,
        as dictionary: key=allele name, value=set of positions

    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __str__(self):
        s = '# Wildtype to be deimmunized\n'
        s += '>' + self.name + '\n' + self.sequence
        s += '\nLength: %i' % len(self.sequence)
        s += '\nOffsets: (uniprot=%i, internal=%i)' % (self.uniprot_offset, self.internal_offset)
        return s

    @classmethod
    def create(cls, alignment, ev_couplings):
        """Construct wildtype from the first sequence in the given alignment.

        Parameters
        ----------
        alignment : `MultipleSeqAlignment`
            Alignment, used for plmc calculation in focus mode
            (with upper and lower letters).
        ev_couplings : `EVcouplings`
            e_ij binaries read in, generated from `alignment`.

        Returns
        -------
        wildtype : `Wildtype`
            The wildtype initialized with its name and sequence.

        """
        wildtype = cls(alignment[0].id, str(alignment[0].seq))
        wildtype.calc_offsets(ev_couplings)
        print '%s\n' % wildtype
        return wildtype

    def calc_offsets(self, ev_couplings):
        """Calculate uniprot and internal offset of the wildtype sequence.

        Parameters
        ----------
        ev_couplings : `EVcouplings`
            e_ij binaries read in, generated from `alignment`.

        """
        # uniprot offset
        self.index_list = ev_couplings.index_list
        self.uniprot_offset = self.index_list[0]
        # internal offset
        self.internal_offset = 0
        for i in xrange(len(self.sequence)):
            if self.sequence[i].isupper():
                self.internal_offset = i
                break
        self.sequence = self.sequence.upper()

    def predict_epitope_positions(self, allele_coll, epitope_length):
        """Predict epitope positions within the wildtype sequence.

        Parameters
        ----------
        allele_coll : `AlleleCollection`
            Collection of alleles for which epitopes are to be predictet.
        epitope_length : int
            Length of an epitope.

        """
        pssms = {(allele.name, k[0], k[1]): v
                 for allele in allele_coll.alleles
                 for k, v in allele.pssm.iteritems()}

        # get starting positions of epitopes (with specified alleles)
        pool = Pool(cpu_count())
        self.epitope_start_pos_alleles = pool.map(
            _predictstart,
            [(_tepitope_pred, (i, self.sequence[i: i + epitope_length], pssms,
                               allele.name, allele.pssm_thresh))
             for i in xrange(len(self.sequence) - (epitope_length - 1))
             for allele in allele_coll.alleles])
        self.epitope_start_pos_alleles = filter(lambda x: x[1] != -1, self.epitope_start_pos_alleles)
        self.num_epitopes = len(self.epitope_start_pos_alleles)

        self.epitope_pos_alleles = defaultdict(set)
        for allele_name, pos, _ in self.epitope_start_pos_alleles:
            for j in range(pos, pos + 9):
                self.epitope_pos_alleles[allele_name].add(j)

        # epitope_pos_ext is the extended epitope set with buffer of 8 aa at the beginning of each epitope
        epitope_start_pos = set([pos for _, pos, _ in self.epitope_start_pos_alleles])
        self.epitope_pos, epitope_pos_ext = [], []
        for i in epitope_start_pos:
            self.epitope_pos.extend(range(i, i + epitope_length))
            if i >= (epitope_length - 2):
                end = i + epitope_length\
                    if i + epitope_length < (len(self.sequence) - 8)\
                    else len(self.sequence) - 8
                epitope_pos_ext.extend(range(i - (epitope_length - 2), end))
            else:
                epitope_pos_ext.extend(range(i + epitope_length))
        self.epitope_pos.sort(), epitope_pos_ext.sort()
        self.epitope_pos, epitope_pos_ext = set(self.epitope_pos), set(epitope_pos_ext)

        # map EVFold uniprot index to internal index
        epitope_pos_mapped = set(self.map_to_seq_idx(pos) for pos in self.index_list)
        if len(epitope_pos_mapped & self.epitope_pos) == 0:
            sys.exit('No epitopes in area considered during Eij calculation! '
                     'Apply normal epitope reduction.')
        self.epitope_pos = epitope_pos_mapped

    def print_epitope_prediction(self, excluded_pos, ignored_pos):
        """For each allele, print epitope prediction of the wildtype sequence.

        Parameters
        ----------
        excluded_pos : list of int
            List of positions not allowed to be mutated.
        ignored_pos : list of int
            List of positions that are fully ignored.

        """
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

    def map_to_seq_idx(self, uniprot_index):
        """Map uniprot index to internal index.

        Parameters
        ----------
        uniprot_index : int
            The uniprot index to be mapped.

        Returns
        -------
        int
            The corresponding internal index.
        """
        return self.internal_offset + (uniprot_index - self.uniprot_offset)

    def map_to_uniprot(self, position):
        """Map position to uniprot index.

        Parameters
        ----------
        position : int
            The sequence position to be mapped.

        Returns
        -------
        int
            The corresponding uniprot index.
        """
        return self.uniprot_offset + (position - self.internal_offset)

    def map_to_seq(self, seq):
        """Map wildtype sequence to another given sequence.

        Parameters
        ----------
        seq : str
            The sequence to be mapped.

        Returns
        -------
        mapping : dict (key=int, value=int)
            Mapping from the wildtype sequence to the given sequence,
            as dictionary: key=residue number of the wildtype,
            value=residue number of the given sequence.
        """
        aln_wt, aln_seq, _, _, _ = pairwise2.align.globalds(
            self.sequence, seq.upper(), blosum62, -11, -1
        )[0]
        mapping = defaultdict(lambda: -1)
        wt_res_num, seq_res_num = 0, 0
        for wt_res, seq_res in zip(aln_wt, aln_seq):
            if wt_res != '-' and seq_res != '-':
                mapping[wt_res_num] = seq_res_num
                wt_res_num += 1
                seq_res_num += 1
            elif wt_res == '-':
                seq_res_num += 1
            elif seq_res == '-':
                wt_res_num += 1
        return mapping

    def to_param_N(self):
        """Return length of the wildtype sequence as string."""
        return str(len(self.sequence))

    def to_set_E(self):
        """Return predicted epitope positions as string."""
        return ' '.join(str(i + 1) for i in self.epitope_pos)
