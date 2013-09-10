import os.path
from collections import namedtuple

def save_trajectory_to_pdb(filename, trajectory, universe, bond_energy):
    with open(filename, 'w') as f:
        for frame_num, frame in enumerate(trajectory):
            f.write("HEADER Coordinates at Frame %d\n" % frame_num)
            f.write("%s" % universe_to_str(universe, frame.coords))
            if frame_num == 0:
                f.write("%s" % connect_records_to_str(bond_energy))
            f.write("END\n")

def universe_to_str(universe, coords):
    my_str = ""
    for i, atom in enumerate(universe):
        x, y, z = (coords[i,0], coords[i,1], coords[i,2])
        my_str += "ATOM  %5d %4s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%6.1f%6.2f%6.2f \n" % \
        ( atom.serial_num, atom.atom_name, atom.res_name, atom.chain_id,
          atom.res_num, x, y, z, atom.charge, 0.00, atom.mass, atom.radius,
          atom.charge )
    return my_str

def connect_records_to_str(bond_energy):
    my_str = ""
    for bond in bond_energy:
        my_str += "CONNECT %d %d\n" % (bond.serial_num_1, bond.serial_num_2)
    return my_str


TrajectoryFrame = namedtuple("TrajectoryFrame", ['coords'])

class Trajectory(object):
    """docstring for Trajectory"""
    def __init__(self):
        super(Trajectory, self).__init__()
        self.frames = []

    def __len__(self):
        return len(self.frames)

    def __iter__(self):
        for frame in self.frames:
            yield frame

    def add_frame(self, coords):
        trajectory_frame = TrajectoryFrame(coords)
        self.frames.append(trajectory_frame)
