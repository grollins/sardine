from numpy import dot, argsort
from scipy.linalg import eig
from numdifftools import Hessian

def compute_hessian(energy_func, X):
    """docstring for compute_hessian"""
    h = Hessian(energy_func)
    return h(X.flatten())

def compute_force_constant_matrix(H, M):
    """docstring for compute_force_constant_matrix"""
    return dot( dot(M, H), M )

def compute_normal_modes(F, discard_trans_and_rot=True):
    eig_vals, eig_vecs = eig(F)
    eig_vals = eig_vals.real
    inds = argsort(eig_vals)

    if discard_trans_and_rot:
        # option to remove translational and rotational modes
        # there are six such modes for for a system with 3N coords
        inds = inds[6:]

    eig_vals = eig_vals[inds]
    eig_vecs = eig_vecs[:,inds]
    modes = NormalModes(eig_vals, eig_vecs)
    return modes


class NormalModes(object):
    """docstring for NormalModes"""
    def __init__(self, frequencies, amplitudes):
        self.frequencies = frequencies
        self.amplitudes = amplitudes

    def get_mode(self, mode_num):
        return self.frequencies[mode_num], self.amplitudes[:,mode_num]

    def get_frequencies(self):
        return self.frequencies

    def get_amplitudes(self):
        return self.amplitudes
