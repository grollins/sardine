from collections import namedtuple

Bond = namedtuple("Bond", ['serial_num_1', 'serial_num_2', 'force_const', 'r_0'])
Angle = namedtuple("Angle", ['serial_num_1', 'serial_num_2', 'serial_num_3',
                             'force_const', 'theta_0'])


class BondFactory(object):
    """docstring for BondFactory"""
    def __init__(self):
        super(BondFactory, self).__init__()

    def create_bond(self, serial_num_1, serial_num_2, force_const, r_0):
        return Bond(serial_num_1=serial_num_1, serial_num_2=serial_num_2,
                    force_const=force_const, r_0=r_0)


class AngleFactory(object):
    """docstring for AngleFactory"""
    def __init__(self):
        super(AngleFactory, self).__init__()

    def create_angle(self, serial_num_1, serial_num_2, serial_num_3,
                     force_const, theta_0):
        return Angle(serial_num_1=serial_num_1, serial_num_2=serial_num_2,
                     serial_num_3=serial_num_3, force_const=force_const,
                     theta_0=theta_0)
