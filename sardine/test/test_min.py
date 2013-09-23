from unittest import TestCase
from os.path import join
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, AngleEnergyFactory
from ..energy import VDWEnergyFactory
from ..energy import EnergyFunctionFactory, GradientFunctionFactory
from ..minimize import BFGSMinimizer
from ..trajectory import save_trajectory_to_pdb


class TestBondMinimization(TestCase):
    def setUp(self):
        pdb_filename = "sardine/test/test_data/triatomic_stretched.pdb"
        sf_filename = "sardine/test/test_data/triatomic.sf"
        uf = UniverseFactory()
        uf.load_atoms_from_file(pdb_filename)
        universe = uf.create_universe()
        self.universe = universe

        bond_energy_factory = BondEnergyFactory()
        bond_energy_factory.load_bonds_from_file(sf_filename)
        self.bond_energy_factory = bond_energy_factory
        bond_energy_func = bond_energy_factory.create_energy_func(num_atoms=len(universe))
        bond_gradient_func = bond_energy_factory.create_gradient_func(num_atoms=len(universe))
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
        save_trajectory_to_pdb('bond_minimization.pdb', trajectory,
                               self.universe, self.bond_energy_factory)
        print "Wrote bond_minimization.pdb"


class TestAngleMinimization(TestCase):
    def setUp(self):
        pdb_filename = "sardine/test/test_data/triatomic_angle.pdb"
        sf_filename = "sardine/test/test_data/triatomic_angle.sf"
        uf = UniverseFactory()
        uf.load_atoms_from_file(pdb_filename)
        universe = uf.create_universe()
        self.universe = universe

        bond_energy_factory = BondEnergyFactory()
        bond_energy_factory.load_bonds_from_file(sf_filename)
        self.bond_energy_factory = bond_energy_factory
        bond_energy_func = bond_energy_factory.create_energy_func(num_atoms=len(universe))
        bond_gradient_func = bond_energy_factory.create_gradient_func(num_atoms=len(universe))
        self.bond_energy_func = bond_energy_func
        self.bond_gradient_func = bond_gradient_func

        angle_energy_factory = AngleEnergyFactory()
        angle_energy_factory.load_angles_from_file(sf_filename)
        self.angle_energy_factory = angle_energy_factory
        angle_energy_func = angle_energy_factory.create_energy_func()
        angle_gradient_func = angle_energy_factory.create_gradient_func()
        self.angle_energy_func = angle_energy_func
        self.angle_gradient_func = angle_gradient_func

        eff = EnergyFunctionFactory()
        eff.add_energy_term('bonds', bond_energy_func)
        eff.add_energy_term('angles', angle_energy_func)
        self.energy_func = eff.create_energy_func(['bonds', 'angles'],
                                                  num_atoms=len(universe))
        
        gff = GradientFunctionFactory()
        gff.add_gradient_term('bonds', bond_gradient_func)
        gff.add_gradient_term('angles', angle_gradient_func)
        self.gradient_func = gff.create_gradient_func(['bonds', 'angles'],
                                                      num_atoms=len(universe))

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
        save_trajectory_to_pdb('angle_minimization.pdb', trajectory,
                               self.universe, self.bond_energy_factory)
        print "Wrote angle_minimization.pdb"


class TestVDWMinimization(TestCase):
    def setUp(self):
        pdb_filename = "sardine/test/test_data/LJ10.pdb"
        sf_filename = "sardine/test/test_data/LJ10.sf"
        uf = UniverseFactory()
        uf.load_atoms_from_file(pdb_filename)
        universe = uf.create_universe()
        self.universe = universe

        vdw_energy_factory = VDWEnergyFactory()
        vdw_energy_factory.load_vdw_from_file(sf_filename)
        self.vdw_energy_factory = vdw_energy_factory
        vdw_energy_func = vdw_energy_factory.create_energy_func()
        vdw_gradient_func = vdw_energy_factory.create_gradient_func()
        self.vdw_energy_func = vdw_energy_func
        self.vdw_gradient_func = vdw_gradient_func

        eff = EnergyFunctionFactory()
        eff.add_energy_term('vdw', vdw_energy_func)
        self.energy_func = eff.create_energy_func(['vdw'], num_atoms=len(universe))
        
        gff = GradientFunctionFactory()
        gff.add_gradient_term('vdw', vdw_gradient_func)
        self.gradient_func = gff.create_gradient_func(['vdw'], num_atoms=len(universe))

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
        save_trajectory_to_pdb('vdw_minimization.pdb', trajectory,
                               self.universe, None)
        print "Wrote vdw_minimization.pdb"
