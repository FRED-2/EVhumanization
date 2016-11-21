from multiprocessing import Pool, cpu_count
from collections import defaultdict


def predictstart(x):
    return predict(*x)


def predict(func, args):
    return func(*args)


def tepitope_pred(pos, ep, pssm, allele, thr):
    is_binder = sum(pssm[allele, ep[i], i] for i in xrange(len(ep))) >= thr
    return (allele,pos) if is_binder else (allele,-1)


class Wildtype(object):
    """Representing the wildtype to be deimmunized.

    attributes:
        name -- name of the wildtype
        sequence -- sequence of the wildtype
        index_list -- TODO
        uniprot_offset -- TODO
        internal_offset -- TODO
        num_epitopes -- number of predicted epitopes
        epitope_pos -- set of positions belonging to an epitope
        epitope_pos_alleles -- positions belonging to an epitope with specified allele,
                               as dictionary: key = name of the allele, value = set of positions
    """

    def __init__(self, name, sequence, index_list=None, uniprot_offset=None,
                 internal_offset=None, num_epitopes=None, epitope_pos=None, epitope_pos_alleles=None):
        self.name = name
        self.sequence = sequence
        self.index_list = index_list
        self.uniprot_offset = uniprot_offset
        self.internal_offset = internal_offset
        self.num_epitopes = num_epitopes
        self.epitope_pos = epitope_pos
        self.epitope_pos_alleles = epitope_pos_alleles

    @classmethod
    def create(cls, alignment, ev_couplings):
        wildtype = cls(alignment[0].id, str(alignment[0].seq))
        wildtype.calc_offsets(ev_couplings)
        print '%s\n' % wildtype
        return wildtype

    def __str__(self):
        s = '# Wildtype to be deimmunized\n'
        s += '>' + self.name + '\n' + self.sequence
        s += '\nLength: %i' % len(self.sequence)
        s += '\nOffsets: (uniprot=%i, internal=%i)' % (self.uniprot_offset, self.internal_offset)
        return s

    def calc_offsets(self, ev_couplings):
        """Calculates uniprot and internal offset of the wildtype sequence."""
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
        """Predicts epitope positions within the wildtype sequence."""
        pssms = {(allele.name, k[0], k[1]):v for allele in allele_coll.alleles for k, v in allele.pssm.iteritems()}

        # get starting positions of epitopes (with specified alleles)
        pool = Pool(cpu_count())
        epitope_start_pos_alleles = pool.map(predictstart, [(tepitope_pred,(i,self.sequence[i:i+9],pssms,allele.name,allele.pssm_thresh))\
            for i in xrange(len(self.sequence) - (epitope_length-1))\
            for allele in allele_coll.alleles])
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
        epitope_pos_mapped = set(self.map_to_seq_idx(pos) for pos in self.index_list)
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
