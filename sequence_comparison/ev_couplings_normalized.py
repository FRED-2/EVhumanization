"""Extension of the `EVcouplings` class.

`ModifiedEVcouplings` additionally stores modified e_ij pair couplings.
Three modification modes are supported:
    (1) MAX: Each e_ij matrix is normalized by its absolute maximum value.
    (2) ABS_MAX: After taking the absolute value of each e_ij entry,
    normalize each (then absolute) e_ij matrix by its maximum value.
    (3) ABS: Set each e_ij entry to its absolute value.

"""

import sys
import numpy as np
from collections import defaultdict
from copy import deepcopy

from utilities.ev_couplings_v4 import EVcouplings


class ModificationMode():

    """Mode used to modify e_ij pair couplings."""

    MAX = 'MAX: Normalize each e_ij matrix by its maximum value.'
    ABS_MAX = 'ABS_MAX: Normalize each absolute e_ij matrix by its maximum value.'
    ABS = 'ABS: Set each e_ij entry to its absolute value.'

    @classmethod
    def getMode(cls, key):
        try:
            return eval('%s.%s' % (cls.__name__, key.upper()))
        except AttributeError as e:
            exit(e.message)


class ModifiedEVcouplings(EVcouplings):

    """Extends EVcouplings class by modified e_ij pair couplings

    Parameters
    ----------
    filename : str
        Binary eij file containing model parameters from plmc software
    mode : str (ModificationMode.{ABS, ABS_MAX, MAX}), optional (default: `ModificationMode.ABS`)
        Modification mode as constant of class `ModificationMode`.
    alphabet : str, optional (default: None)
        Symbols corresponding to model states (e.g. "-ACGT").
    precision : {'float32', 'float64'}, optional (default: 'float32')
        Sets if input file has single (default) or double precision

    Attributes (additional)
    ----------
    mode : str
        Modification mode.
    modified_e_ij : np.array
        Numpy array of shape (L, L, #symbols, #symbols), where L is the model
        length and #symbols is the size of the alphabet. Contains the
        modified e_ij values.

    """

    def __init__(self, filename, mode=ModificationMode.ABS, alphabet=None,
                 precision='float32'):
        super(ModifiedEVcouplings, self).__init__(filename, alphabet, precision)
        self.mode = mode
        self.modify_eijs()

    def modify_eijs(self):
        """Modify eij matrices in a certain way, specified by the mode."""
        self.modified_e_ij = np.zeros(self.e_ij.shape)
        for i in xrange(self.e_ij.shape[0]):
            for j in xrange(i + 1, self.e_ij.shape[1]):
                self.modify_single_eij_matrix(i, j)

    def modify_single_eij_matrix(self, i, j):
        """Modify a single e_ij matrix."""
        eij_matrix = self.e_ij[i, j]
        abs_eij_matrix = abs(eij_matrix)

        if self.mode == ModificationMode.ABS:
            self.modified_e_ij[i, j] = abs_eij_matrix
            return
        elif self.mode == ModificationMode.ABS_MAX:
            base_eij_matrix = abs_eij_matrix
        elif self.mode == ModificationMode.MAX:
            base_eij_matrix = eij_matrix
        else:
            exit('Modification mode %s is not implemented.' % self.mode)

        max_eij = np.amax(abs_eij_matrix)
        norm = np.vectorize(lambda eij_value: eij_value / float(max_eij))
        self.modified_e_ij = normalize(base_eij_matrix)
        return normalize(base_eij_matrix)
