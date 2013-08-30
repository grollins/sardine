from scipy.spatial.distance import pdist, squareform

def compute_distance_vector(X):
    """docstring for compute_distance_matrix"""
    return pdist(X, 'euclidean')

def compute_distance_matrix_from_vector(D_vec):
    """docstring for compute_distance_matrix"""
    return squareform( D_vec )
