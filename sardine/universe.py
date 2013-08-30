from numpy import array, zeros, zeros_like, sqrt, diagflat
from collections import namedtuple

Atom = namedtuple("Atom", ['x', 'y', 'z', 'mass', 'charge', 'radius',
                           'serial_num', 'res_num', 'atom_name', 'res_name',
                           'chain_id'])


class UniverseFactory(object):
    """docstring for UniverseFactory"""
    def __init__(self):
        self.atoms = []

    def add_atom(self, x, y, z, mass, charge, radius, serial_num=None,
                 res_num='1', atom_name='LJ', res_name='A', chain_id='A'):
        """docstring"""
        serial_num = (len(self.atoms) if serial_num is None else serial_num)
        self.atoms.append( Atom(x, y, z, mass, charge, radius,
                                serial_num, res_num, atom_name, res_name,
                                chain_id) )

    def create_universe(self):
        return Universe(self.atoms)


class Universe(object):
    """docstring for Universe"""
    def __init__(self, atoms):
        self.atoms = atoms
        self.initialize_matrices()

    def __len__(self):
        return len(self.atoms)

    def initialize_matrices(self):
        masses = [a.mass for a in self.atoms]
        coords = [(a.x, a.y, a.z) for a in self.atoms]
        charges = [a.charge for a in self.atoms]
        radii = [a.radius for a in self.atoms]
        self.mass_array = array(masses)
        self.coord_array = array(coords)
        self.charge_array = array(charges)
        self.radius_array = array(radii)
        inv_sqrt_diag = 1./sqrt(self.mass_array)
        self.M = diagflat(inv_sqrt_diag)

    def get_inv_sqrt_mass_matrix(self):
        return self.M

    def get_coords(self):
        return self.coord_array
