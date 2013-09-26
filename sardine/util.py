from numpy import arccos, radians, degrees, pi, dot, isnan, allclose, seterr
from numpy.linalg import norm
from scipy.spatial.distance import pdist, squareform

# numpy complains when arg to arccos is out of bounds, but the
# isnan conditional takes care of that, so we'll supress these warnings
seterr(invalid='ignore')

def compute_distance_vector(X):
    """docstring for compute_distance_matrix"""
    return pdist(X, 'euclidean')

def compute_distance_matrix_from_vector(D_vec):
    """docstring for compute_distance_matrix"""
    return squareform( D_vec )

def compute_angle(u, v):
    """docstring for compute_angle"""
    cos_theta = dot(u, v) / (norm(u) * norm(v))
    theta = arccos(cos_theta)

    if isnan(theta):
        if allclose(u,v):
            theta = 0.0
        else:
            theta = pi
    return theta

def deg2rad(angle):
    """docstring for deg2rad"""
    return radians(angle)

def rad2deg(angle):
    """docstring for rad2deg"""
    return degrees(angle)

def coords_2d_to_1d(X):
    return X.flatten()

def coords_1d_to_2d(X):
    num_atoms = len(X) / 3
    return X.reshape((num_atoms, 3))
