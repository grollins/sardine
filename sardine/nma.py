from numpy import dot
from numdifftools import Hessian
from .util import compute_distance_vector

def compute_hessian(energy_func, X):
    """docstring for compute_hessian"""
    h = Hessian(energy_func)
    D_vec = compute_distance_vector(X)
    return h(D_vec)

def compute_force_constant_matrix(H, M):
    """docstring for compute_force_constant_matrix"""
    return dot( dot(M, H), M )

