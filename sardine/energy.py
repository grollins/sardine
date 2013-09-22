from collections import namedtuple
from numpy import zeros, zeros_like, triu_indices, seterr, dot, pi, cos
from numpy.linalg import norm
from .util import compute_distance_vector, compute_distance_matrix_from_vector,\
                  compute_angle
from .parsers import StructureParser

# numpy complains about division by zero when computing vdw energy because
# the distance matrix has zeros along the diagonal. We only include the terms
# above the diagonal when computing the energy, so we'll suppress these
# warnings.
seterr(divide='ignore')

class BondEnergyFactory(object):
    """docstring for BondEnergyFactory"""
    def __init__(self):
        self.bonds = []
        self.sf_parser = StructureParser()

    def __len__(self):
        return len(self.bonds)

    def __iter__(self):
        for b in self.bonds:
            yield b

    def add_bond(self, bond):
        self.bonds.append( bond )

    def load_bonds_from_file(self, filename):
        if filename.endswith(".sf"):
            for b in self.sf_parser.iter_bonds_in_sf_file(filename):
                self.add_bond(b)
        else:
            print "Expected a .sf file, got %s" % filename
            return

    def create_energy_func(self, num_atoms):
        K = zeros([num_atoms, num_atoms])
        R = zeros([num_atoms, num_atoms])
        for b in self.bonds:
            i = b.serial_num_1 - 1 # convert to 0-indexing
            j = b.serial_num_2 - 1 # convert to 0-indexing
            K[i,j] = b.force_const
            R[i,j] = b.r_0
        def bond_energy_func(X, D):
            E = 0.5 * K * (D - R)**2
            return E.sum()
        return bond_energy_func

    def create_gradient_func(self, num_atoms):
        K = zeros([num_atoms, num_atoms])
        R = zeros([num_atoms, num_atoms])
        for b in self.bonds:
            i = b.serial_num_1 - 1 # convert to 0-indexing
            j = b.serial_num_2 - 1 # convert to 0-indexing
            K[i,j] = b.force_const
            R[i,j] = b.r_0

        def bond_gradient_func(X, D):
            F = -K * (D - R)
            G = zeros_like(X)
            for b in self.bonds:
                i = b.serial_num_1 - 1 # convert to 0-indexing
                j = b.serial_num_2 - 1 # convert to 0-indexing
                dX = (X[j,:] - X[i,:]) / D[j,i]
                G[i,:] += F[i,j] * dX
                G[j,:] += -F[i,j] * dX
            return G
        return bond_gradient_func


class AngleEnergyFactory(object):
    """docstring for AngleEnergyFactory"""
    def __init__(self):
        super(AngleEnergyFactory, self).__init__()
        self.angles = []
        self.sf_parser = StructureParser()

    def __len__(self):
        return len(self.angles)

    def __iter__(self):
        for a in self.angles:
            yield a

    def add_angle(self, angle):
        self.angles.append( angle )

    def load_angles_from_file(self, filename):
        if filename.endswith(".sf"):
            for a in self.sf_parser.iter_angles_in_sf_file(filename):
                self.add_angle(a)
        else:
            print "Expected a .sf file, got %s" % filename
            return

    def create_energy_func(self):
        def angle_energy_func(X, D):
            E = 0.0
            for a in self.angles:
                i = a.serial_num_1 - 1 # convert to 0-indexing
                j = a.serial_num_2 - 1 # convert to 0-indexing
                k = a.serial_num_3 - 1 # convert to 0-indexing
                vec12 = X[i,:] - X[j,:]
                vec32 = X[k,:] - X[j,:]
                theta = compute_angle(vec12, vec32)
                E += 0.5 * a.force_const * (theta - a.theta_0)**2
            return E.sum()
        return angle_energy_func

    def create_gradient_func(self):
        def angle_gradient_func(X, D):
            G = zeros_like(X)
            for a in self.angles:
                i = a.serial_num_1 - 1 # convert to 0-indexing
                j = a.serial_num_2 - 1 # convert to 0-indexing
                k = a.serial_num_3 - 1 # convert to 0-indexing
                vec12 = X[i,:] - X[j,:]
                vec32 = X[k,:] - X[j,:]
                theta = compute_angle(vec12, vec32)
                total_G = a.force_const * (theta - a.theta_0)

                length_vec12 = norm(vec12)
                length_vec32 = norm(vec32)
                vec12 /= length_vec12
                vec32 /= length_vec32
                cos_theta = cos(theta)
                G1 = total_G * (vec32 - cos_theta * vec12) / length_vec12
                G2 = total_G * (vec12 - cos_theta * vec32) / length_vec32
                G_atom1 = -G1
                G_atom2 = G1 + G2
                G_atom3 = -G2

                G = zeros_like(X)
                G[i,:] = G_atom1
                G[j,:] = G_atom2
                G[k,:] = G_atom3
            return G
        return angle_gradient_func


class VDWEnergyFactory(object):
    """docstring for VDWEnergyFactory"""
    def __init__(self):
        super(VDWEnergyFactory, self).__init__()
        self.well_distance = 2.6 # angstroms
        self.well_depth = 0.1 # kcal/mol
        self.sf_parser = StructureParser()

    def set_well_distance(self, well_distance):
        self.well_distance = well_distance

    def set_well_depth(self, well_depth):
        self.well_depth = well_depth

    def load_vdw_from_file(self, filename):
        if filename.endswith(".sf"):
            # take parameters from first VDW line encountered in file
            distance, depth = self.sf_parser.get_first_vdw_in_sf_file(filename)
            self.set_well_distance(distance)
            self.set_well_depth(depth)
        else:
            print "Expected a .sf file, got %s" % filename
            return

    def create_energy_func(self):
        well_distance = self.well_distance
        well_depth = self.well_depth
        def vdw_energy_func(X, D):
            separation_ratio = well_distance / D
            sr6 = separation_ratio ** 6
            sr12 = sr6 * sr6
            E = well_depth * (sr12 - (2 * sr6))
            # only sum terms above diagonal to avoid double counting
            inds = triu_indices( E.shape[0], 1 )
            return E[inds].sum()
        return vdw_energy_func


class EnergyFunctionFactory(object):
    """docstring for EnergyFunctionFactory"""
    def __init__(self):
        self.energy_terms = {}

    def add_energy_term(self, term_name, term_func):
        self.energy_terms[term_name] = term_func

    def create_energy_func(self, term_names, num_atoms):
        def energy_func(X_vec):
            X = X_vec.reshape((num_atoms, 3))
            D_vec = compute_distance_vector(X)
            D = compute_distance_matrix_from_vector(D_vec)
            E = 0.
            for t in term_names:
                E += self.energy_terms[t](X, D)
            return E
        return energy_func


class GradientFunctionFactory(object):
    """docstring for GradientFunctionFactory"""
    def __init__(self):
        self.gradient_terms = {}

    def add_gradient_term(self, term_name, term_func):
        self.gradient_terms[term_name] = term_func

    def create_gradient_func(self, term_names, num_atoms):
        def gradient_func(X_vec):
            X = X_vec.reshape((num_atoms, 3))
            D_vec = compute_distance_vector(X)
            D = compute_distance_matrix_from_vector(D_vec)
            G = zeros_like(X)
            for t in term_names:
                G += self.gradient_terms[t](X,D)
            return G.flatten()
        return gradient_func
