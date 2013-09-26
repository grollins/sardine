"""
Computes normal modes for a triatomic molecule and compares them to
analytical results presented in
Molecular Modeling: Principles and Applications by A.R. Leach, pg. 275 (2nd ed.)
"""
from os.path import exists
from unittest import TestCase
from numpy import array, sqrt, dot, diagflat, allclose, zeros, repeat
from ..universe import UniverseFactory
from ..energy import BondEnergyFactory, EnergyFunctionFactory
from ..nma import compute_hessian, compute_force_constant_matrix,\
                  compute_normal_modes
from ..const import CM_CONVERSION_FACTOR

PDB_FILENAME = "sardine/test/test_data/triatomic.pdb"
SF_FILENAME = "sardine/test/test_data/triatomic.sf"


class TestTriatomicNMA(TestCase):
    def setUp(self):
        m1 = 1.0
        m2 = 1.0
        m3 = 1.0
        k = 0.1
        r_0 = 1.0

        uf = UniverseFactory()
        uf.load_atoms_from_file(PDB_FILENAME)
        universe = uf.create_universe()
        self.universe = universe

        bond_energy = BondEnergyFactory()
        bond_energy.load_bonds_from_file(SF_FILENAME)
        bond_energy_func = bond_energy.create_energy_func(num_atoms=len(universe))
        self.bond_energy_func = bond_energy_func
        eff = EnergyFunctionFactory()
        eff.add_energy_term('bonds', bond_energy_func)
        self.energy_func = eff.create_energy_func(['bonds'], num_atoms=len(universe))

        H = zeros((9,9))
        H[0,0] = k
        H[0,3] = -k
        H[3,0] = -k
        H[3,3] = 2*k
        H[3,6] = -k
        H[6,0] = 0
        H[6,3] = -k
        H[6,6] = k
        expected_H = H

        M_diag = array((m1, m2, m3))
        num_spatial_dimensions = 3
        inv_sqrt_diag = 1./sqrt(repeat(M_diag, num_spatial_dimensions))
        expected_M = diagflat(inv_sqrt_diag)

        expected_F = dot( dot(expected_M, expected_H), expected_M )
        self.expected_mode_freqs = array([0.0, k/m1, k*(m2+2*m1)/(m1*m2)])
        self.expected_mode_freqs = sqrt(self.expected_mode_freqs) * CM_CONVERSION_FACTOR
        self.expected_H = expected_H
        self.expected_M = expected_M
        self.expected_F = expected_F

    def test_computes_correct_normal_modes(self):
        M = self.universe.get_inv_sqrt_mass_matrix()
        self.assertTrue( allclose(M, self.expected_M),
                         msg="\n%s\n%s" % (M, self.expected_M) )

        X = self.universe.get_coords()
        H = compute_hessian(self.energy_func, X)
        self.assertTrue( allclose(H, self.expected_H),
                         msg="\n%s\n%s" % (H, self.expected_H) )

        F = compute_force_constant_matrix(H, M)
        self.assertTrue( allclose(F, self.expected_F),
                         msg="\n%s\n%s" % (F, self.expected_F) )

        normal_modes = compute_normal_modes(F, discard_trans_and_rot=True)
        mode_freqs = normal_modes.get_frequencies()
        print mode_freqs
        self.assertTrue( allclose(mode_freqs, self.expected_mode_freqs, atol=1e-04),
                         msg="\n%s\n%s" % (mode_freqs, self.expected_mode_freqs) )
