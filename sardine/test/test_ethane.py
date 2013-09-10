"""
Computes normal modes for ethane.
"""

from unittest import TestCase
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, AngleEnergyFactory, EnergyFunctionFactory
from ..nma import compute_hessian, compute_force_constant_matrix,\
                  compute_normal_modes

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

        angle_energy = AngleEnergyFactory()
        angle_energy.load_angles_from_file(SF_FILENAME)
        self.assertEqual( len(angle_energy), 12 )
        angle_energy_func = angle_energy.create_energy_func()

        eff = EnergyFunctionFactory()
        eff.add_energy_term('bonds', bond_energy_func)
        eff.add_energy_term('angles', angle_energy_func)
        energy_func = eff.create_energy_func(
                        ['bonds', 'angles'], num_atoms=len(universe))
        self.universe = universe
        self.energy_func = energy_func

    def test_computes_normal_modes(self):
        M = self.universe.get_inv_sqrt_mass_matrix()
        X = self.universe.get_coords()
        H = compute_hessian(self.energy_func, X)
        F = compute_force_constant_matrix(H, M)
        normal_modes = compute_normal_modes(F, discard_trans_and_rot=False)
        # print normal_modes
