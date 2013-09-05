from collections import namedtuple

Atom = namedtuple("Atom", ['x', 'y', 'z', 'mass', 'charge', 'radius',
                           'serial_num', 'res_num', 'atom_name', 'res_name',
                           'chain_id'])

class AtomFactory(object):
    def __init__(self):
        super(AtomFactory, self).__init__()

    def create_atom(self, x, y, z, mass, charge, radius, serial_num=None,
                    res_num='1', atom_name='LJ', res_name='A', chain_id='A'):
        return Atom(x, y, z, mass, charge, radius,
                    serial_num, res_num, atom_name, res_name, chain_id)
