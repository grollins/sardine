from unittest import TestCase
from os.path import join
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, EnergyFunctionFactory
from ..energy import GradientFunctionFactory
from ..minimize import BFGSMinimizer
from ..trajectory import save_trajectory_to_pdb


PDB_FILENAME = "sardine/test/test_data/triatomic_stretched.pdb"
SF_FILENAME = "sardine/test/test_data/triatomic.sf"


class TestBondMinimization(TestCase):
    def setUp(self):
        uf = UniverseFactory()
        uf.load_atoms_from_file(PDB_FILENAME)
        universe = uf.create_universe()
        self.universe = universe

        bond_energy = BondEnergyFactory()
        bond_energy.load_bonds_from_file(SF_FILENAME)
        self.bond_energy = bond_energy
        bond_energy_func = bond_energy.create_energy_func(num_atoms=len(universe))
        bond_gradient_func = bond_energy.create_gradient_func(num_atoms=len(universe))
        self.bond_energy_func = bond_energy_func
        self.bond_gradient_func = bond_gradient_func

        eff = EnergyFunctionFactory()
        eff.add_energy_term('bonds', bond_energy_func)
        self.energy_func = eff.create_energy_func(['bonds'], num_atoms=len(universe))
        
        gff = GradientFunctionFactory()
        gff.add_gradient_term('bonds', bond_gradient_func)
        self.gradient_func = gff.create_gradient_func(['bonds'], num_atoms=len(universe))

    def test_minimizes_to_expected_coordinates(self):
        minimizer = BFGSMinimizer()
        X = self.universe.get_coords().flatten()
        energy_initial = self.energy_func(X)
        X_min, energy_min = minimizer.run_minimization(
                                self.energy_func, self.gradient_func, X,
                                num_atoms=len(self.universe),
                                save_trajectory=True, noisy=True)
        print energy_initial, energy_min
        print X_min

        trajectory = minimizer.get_trajectory()
        save_trajectory_to_pdb('triatomic_minimization.pdb', trajectory,
                               self.universe, self.bond_energy)
        print "Wrote triatomic_minimization.pdb"
