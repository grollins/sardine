from collections import namedtuple
from numpy import zeros
from .util import compute_distance_matrix_from_vector

Bond = namedtuple("Bond", ['serial_num_1', 'serial_num_2', 'force_const', 'r_0'])

class BondEnergyFactory(object):
    """docstring for BondEnergy"""
    def __init__(self):
        self.bonds = []

    def load_structure_file(self, structure_file):
        with open(structure_file, 'r') as f:
            for line in f:
                word_list = line.split()
                if word_list[0] == 'BOND':
                    serial_num_1 = int(word_list[1])
                    serial_num_2 = int(word_list[2])
                    force_const = float(word_list[3])
                    r_0 = float(word_list[4])
                    self.add_term(serial_num_1, serial_num_2, force_const, r_0)

    def add_term(self, serial_num_1, serial_num_2, force_const, r_0):
        self.bonds.append( Bond(serial_num_1, serial_num_2, force_const, r_0) )

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

    def create_energy_func(self, term_names):
        def energy_func(D_vec):
            D = compute_distance_matrix_from_vector(D_vec)
            E = 0.
            for t in term_names:
                E += self.energy_terms[t](D)
            return E
        return energy_func
