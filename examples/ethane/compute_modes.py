"""
Computes normal modes for ethane.
"""

from os import mkdir
from os.path import exists, join
from sardine.universe import UniverseFactory
from sardine.energy import BondEnergyFactory, AngleEnergyFactory, VDWEnergyFactory
from sardine.energy import EnergyFunctionFactory, GradientFunctionFactory
from sardine.nma import compute_hessian, compute_force_constant_matrix,\
                        compute_normal_modes, generate_mode_trajectory
from sardine.trajectory import save_trajectory_to_pdb
from sardine.minimize import BFGSMinimizer


PDB_FILENAME = "C2H6_ideal_trans_min_final.pdb"
SF_FILENAME = "C2H6.sf"
OUTPUT_DIR = "modes"


def main():
    if not exists(OUTPUT_DIR):
        mkdir(OUTPUT_DIR)

    uf = UniverseFactory()
    uf.load_atoms_from_file(PDB_FILENAME)
    universe = uf.create_universe()

    bond_energy_factory = BondEnergyFactory()
    bond_energy_factory.load_bonds_from_file(SF_FILENAME)
    bond_energy_func = bond_energy_factory.create_energy_func(num_atoms=len(universe))
    bond_gradient_func = bond_energy_factory.create_gradient_func(num_atoms=len(universe))

    angle_energy_factory = AngleEnergyFactory()
    angle_energy_factory.load_angles_from_file(SF_FILENAME)
    angle_energy_func = angle_energy_factory.create_energy_func()
    angle_gradient_func = angle_energy_factory.create_gradient_func()

    vdw_energy_factory = VDWEnergyFactory()
    vdw_energy_factory.load_vdw_from_file(SF_FILENAME)
    vdw_energy_func = vdw_energy_factory.create_energy_func()
    vdw_gradient_func = vdw_energy_factory.create_gradient_func()

    eff = EnergyFunctionFactory()
    eff.add_energy_term('bonds', bond_energy_func)
    eff.add_energy_term('angles', angle_energy_func)
    eff.add_energy_term('vdw', vdw_energy_func)
    energy_func = eff.create_energy_func(
                    ['bonds', 'angles', 'vdw'], num_atoms=len(universe))

    gff = GradientFunctionFactory()
    gff.add_gradient_term('bonds', bond_gradient_func)
    gff.add_gradient_term('angles', angle_gradient_func)
    gff.add_gradient_term('vdw', vdw_gradient_func)
    gradient_func = gff.create_gradient_func(
                        ['bonds', 'angles', 'vdw'], num_atoms=len(universe))

    # ======================
    # = Minimize structure =
    # ======================
    minimizer = BFGSMinimizer(maxiter=200)
    X = universe.get_coords().flatten()
    energy_initial = energy_func(X)
    X_min, energy_min = minimizer.run_minimization(
                            energy_func, gradient_func, X,
                            num_atoms=len(universe),
                            save_trajectory=True, noisy=True)
    print energy_initial, energy_min
    trajectory = minimizer.get_trajectory()
    save_trajectory_to_pdb('minimization.pdb', trajectory, universe,
                           bond_energy_factory)
    print "Wrote minimization.pdb"

    # ========================
    # = Compute normal modes =
    # ========================
    M = universe.get_inv_sqrt_mass_matrix()
    # X = universe.get_coords()
    X = X_min # use minimized coordinates
    H = compute_hessian(energy_func, X)
    F = compute_force_constant_matrix(H, M)
    normal_modes = compute_normal_modes(F, discard_trans_and_rot=True)
    mode_freqs = normal_modes.get_frequencies()

    with open(join(OUTPUT_DIR, 'eigen_values.txt'), 'w') as f:
        f.write("%s" % normal_modes.freq_to_str())

    for i in xrange(len(mode_freqs)):
        mode_trajectory = generate_mode_trajectory(universe, normal_modes,
                                                   mode_number=i)
        save_trajectory_to_pdb(
            join(OUTPUT_DIR, 'ethane_mode%02d.pdb') % (i+1),
            mode_trajectory, universe, bond_energy_factory)

if __name__ == '__main__':
    main()