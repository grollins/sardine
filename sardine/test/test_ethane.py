"""
Computes normal modes for ethane.
"""

from unittest import TestCase
from os.path import exists
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, AngleEnergyFactory, VDWEnergyFactory,\
                     EnergyFunctionFactory
from ..nma import compute_hessian, compute_force_constant_matrix,\
                  compute_normal_modes, generate_mode_trajectory
from ..trajectory import save_trajectory_to_pdb

PDB_FILENAME = "sardine/test/test_data/C2H6_ideal_trans_min_final.pdb"
SF_FILENAME = "sardine/test/test_data/C2H6.sf"

class TestEthaneNMA(TestCase):
    def setUp(self):
        uf = UniverseFactory()
        uf.load_atoms_from_file(PDB_FILENAME)
        universe = uf.create_universe()

        bond_energy = BondEnergyFactory()
        bond_energy.load_bonds_from_file(SF_FILENAME)
        self.assertEqual( len(bond_energy), 7 )
        bond_energy_func = bond_energy.create_energy_func(num_atoms=len(universe))
        self.bond_energy = bond_energy

        angle_energy = AngleEnergyFactory()
        angle_energy.load_angles_from_file(SF_FILENAME)
        self.assertEqual( len(angle_energy), 12 )
        angle_energy_func = angle_energy.create_energy_func()

        vdw_energy = VDWEnergyFactory()
        vdw_energy.load_vdw_from_file(SF_FILENAME)
        vdw_energy_func = vdw_energy.create_energy_func()

        eff = EnergyFunctionFactory()
        eff.add_energy_term('bonds', bond_energy_func)
        eff.add_energy_term('angles', angle_energy_func)
        eff.add_energy_term('vdw', vdw_energy_func)
        energy_func = eff.create_energy_func(
                        ['bonds', 'angles', 'vdw'], num_atoms=len(universe))
        self.universe = universe
        self.energy_func = energy_func

    def test_computes_normal_modes(self):
        M = self.universe.get_inv_sqrt_mass_matrix()
        X = self.universe.get_coords()
        H = compute_hessian(self.energy_func, X)
        F = compute_force_constant_matrix(H, M)
        normal_modes = compute_normal_modes(F, discard_trans_and_rot=False)
        mode_freqs = normal_modes.get_frequencies()
        print mode_freqs
        for i in xrange(len(mode_freqs)):
            mode_trajectory = generate_mode_trajectory(self.universe, normal_modes,
                                                       mode_number=i)
            print len(mode_trajectory), "frames"
            save_trajectory_to_pdb('ethane_traj_mode%02d.pdb' % (i+1),
                                   mode_trajectory,
                                   self.universe, self.bond_energy)
            self.assertTrue( exists('ethane_traj_mode%02d.pdb' % (i+1)) )
