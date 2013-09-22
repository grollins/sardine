from .atom import AtomFactory
from .bonded_terms import BondFactory, AngleFactory
from .util import deg2rad

class PdbParser(object):
    """
    Expected format
    `ATOM serial_num atom_name res_name chain_id res_num x y z unused unused mass radius charge`
    Example
    `ATOM 1 C A A 1 -0.002 0.000 0.000 0.00 1.00 12.0 1.75 0.00`
    """
    def __init__(self):
        super(PdbParser, self).__init__()
        self.atom_factory = AtomFactory()

    def iter_atoms_in_pdb_file(self, filename):
        with open(filename, 'r') as f:
            lines = [line.split() for line in f.readlines()]
            for L in lines:
                if L[0] == 'ATOM':
                    serial_num = int(L[1])
                    atom_name = L[2]
                    res_name = L[3]
                    chain_id = L[4]
                    res_num = int(L[5])
                    x = float(L[6])
                    y = float(L[7])
                    z = float(L[8])
                    unused = L[9]
                    unused = L[10]
                    mass = float(L[11])
                    radius = float(L[12])
                    charge = float(L[13])
                    this_atom = self.atom_factory.create_atom(
                                    x=x, y=y, z=z, mass=mass, charge=charge,
                                    radius=radius, serial_num=serial_num,
                                    res_num=res_num, atom_name=atom_name,
                                    res_name=res_name, chain_id=chain_id)
                    yield this_atom
                else:
                    continue


class StructureParser(object):
    """docstring for StructureParser"""
    def __init__(self):
        super(StructureParser, self).__init__()
        self.bond_factory = BondFactory()
        self.angle_factory = AngleFactory()
        # self.torsion_factory = TorsionFactory()

    def iter_bonds_in_sf_file(self, filename):
        with open(filename, 'r') as f:
            lines = [line.split() for line in f.readlines()]
            for L in lines:
                if L[0] == 'BOND':
                    serial_num_1 = int(L[1])
                    serial_num_2 = int(L[2])
                    force_const = float(L[3])
                    r_0 = float(L[4])
                    this_bond = self.bond_factory.create_bond(
                                    serial_num_1=serial_num_1,
                                    serial_num_2=serial_num_2,
                                    force_const=force_const, r_0=r_0)
                    yield this_bond
                else:
                    continue

    def iter_angles_in_sf_file(self, filename):
        with open(filename, 'r') as f:
            lines = [line.split() for line in f.readlines()]
            for L in lines:
                if L[0] == 'ANGLE':
                    serial_num_1 = int(L[1])
                    serial_num_2 = int(L[2])
                    serial_num_3 = int(L[3])
                    force_const = float(L[4])
                    theta_0 = deg2rad( float(L[5]) )
                    this_angle = self.angle_factory.create_angle(
                                    serial_num_1=serial_num_1,
                                    serial_num_2=serial_num_2,
                                    serial_num_3=serial_num_3,
                                    force_const=force_const,
                                    theta_0=theta_0)
                    yield this_angle
                else:
                    continue

    def get_first_vdw_in_sf_file(self, filename):
        with open(filename, 'r') as f:
            lines = [line.split() for line in f.readlines()]
            for L in lines:
                if L[0] == 'VDW':
                    well_distance = float(L[1])
                    well_depth = float(L[2])
                    return (well_distance, well_depth)
                else:
                    continue
