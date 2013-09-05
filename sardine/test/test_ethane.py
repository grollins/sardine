"""
Computes normal modes for ethane.
"""

from unittest import TestCase
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, EnergyFunctionFactory
from ..nma import compute_hessian, compute_force_constant_matrix,\
                  compute_normal_modes

PDB_FILENAME = "sardine/test/test_data/C2H6_ideal_trans_min_final.pdb"
SF_FILENAME = "sardine/test/test_data/C2H6.sf"

class TestEthaneNMA(TestCase):
    def setUp(self):
        uf = UniverseFactory()
        uf.load_atoms_from_file(PDB_FILENAME)
        universe = uf.create_universe()

        bef = BondEnergyFactory()
        bef.load_bonds_from_file(SF_FILENAME)
        bond_energy_func = bef.create_func(num_atoms=len(universe))
        eff = EnergyFunctionFactory()
        eff.add_energy_term('bonds', bond_energy_func)
        energy_func = eff.create_energy_func(
                        ['bonds'], num_atoms=len(universe))
        self.universe = universe
        self.energy_func = energy_func

    def test_computes_normal_modes(self):
        M = self.universe.get_inv_sqrt_mass_matrix()
        X = self.universe.get_coords()
        H = compute_hessian(self.energy_func, X)
        F = compute_force_constant_matrix(H, M)
        normal_modes = compute_normal_modes(F, discard_trans_and_rot=False)
        print normal_modes
