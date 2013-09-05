from numpy import array, zeros, zeros_like, sqrt, diagflat, repeat
from .parsers import PdbParser

class UniverseFactory(object):
    """docstring for UniverseFactory"""
    def __init__(self):
        self.atoms = []
        self.pdb_parser = PdbParser()

    def add_atom(self, atom):
        """docstring"""
        if atom.serial_num is None:
            atom.serial_num = len(self.atoms)
        self.atoms.append(atom)

    def create_universe(self):
        return Universe(self.atoms)

    def load_atoms_from_file(self, filename):
        if filename.endswith(".pdb"):
            for a in self.pdb_parser.iter_atoms_in_pdb_file(filename):
                self.add_atom(a)
        else:
            print "Expected a .pdb file, got %s" % filename
            return

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
        inv_sqrt_diag = 1./sqrt( repeat(self.mass_array, 3) )
        self.M = diagflat(inv_sqrt_diag)

    def get_inv_sqrt_mass_matrix(self):
        return self.M

    def get_coords(self):
        return self.coord_array
