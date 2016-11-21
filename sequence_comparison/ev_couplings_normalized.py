import sys
import numpy as np
from collections import defaultdict
from copy import deepcopy

from ev_couplings_v4 import EVcouplings

class NormalizedEVcouplings(EVcouplings):
    """Extends EVcouplings class by normalized eij matrices

    Additional attributes:
        mode -- normalization mode
        normalized_e_ij -- normalized eij matrices
    """

    # normalize each eij matrix by its maximum value
    MAX = 'max'

    # normalize each absolute eij matrix by its maximum value
    ABS_MAX = 'abs max'

    def __init__(self, filename, mode=ABS_MAX,
                 alphabet=None, precision="float32"):
        super(NormalizedEVcouplings, self).__init__(filename, alphabet, precision)
        self.mode = mode
        self.normalize_eijs(mode)

    def normalize_eijs(self, mode=ABS_MAX):
        """Normalize eij matrices in a certain way, specified by the mode."""
        self.normalized_e_ij = np.zeros(self.e_ij.shape)
        for i in xrange(self.e_ij.shape[0]):
            for j in xrange(i+1, self.e_ij.shape[1]):
                self.normalized_e_ij[i, j]\
                    = NormalizedEVcouplings.normalize_eij(self.e_ij[i, j], mode)

    @staticmethod
    def normalize_eij(eij, mode):
        """Normalize a single eij matrix in a certain way, specified by the mode."""
        if mode in [NormalizedEVcouplings.MAX, NormalizedEVcouplings.ABS_MAX]:
            abs_eij = abs(eij)
            norm_factor = np.amax(abs_eij)
        if mode == NormalizedEVcouplings.MAX:
            pass
        elif mode == NormalizedEVcouplings.ABS_MAX:
            eij = abs_eij
        else:
            print >> sys.stderr, 'Mode %s is not implemented.' % mode
            exit(1)
        norm = np.vectorize(lambda x: x / float(norm_factor))
        return norm(eij)
