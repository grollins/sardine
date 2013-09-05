from collections import namedtuple
from numpy import zeros
from .util import compute_distance_vector, compute_distance_matrix_from_vector
from .parsers import StructureParser

class BondEnergyFactory(object):
    """docstring for BondEnergy"""
    def __init__(self):
        self.bonds = []
        self.sf_parser = StructureParser()

    def add_bond(self, bond):
        self.bonds.append( bond )

    def load_bonds_from_file(self, filename):
        if filename.endswith(".sf"):
            for b in self.sf_parser.iter_bonds_in_sf_file(filename):
                self.add_bond(b)
        else:
            print "Expected a .sf file, got %s" % filename
            return

    def create_func(self, num_atoms):
        K = zeros([num_atoms, num_atoms])
        R = zeros([num_atoms, num_atoms])
        for b in self.bonds:
            i = b.serial_num_1 - 1 # convert to 0-indexing
            j = b.serial_num_2 - 1 # convert to 0-indexing
            K[i,j] = b.force_const
            R[i,j] = b.r_0

        def bond_energy_func(D):
            E = 0.5 * K * (D - R)**2
            return E.sum()

        return bond_energy_func


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
                E += self.energy_terms[t](D)
            return E
        return energy_func
