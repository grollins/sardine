from numpy import dot, argsort, linspace, concatenate
from scipy.linalg import eig
from numdifftools import Hessian
from .trajectory import Trajectory

def compute_hessian(energy_func, X):
    """docstring for compute_hessian"""
    h = Hessian(energy_func)
    return h(X.flatten())

def compute_force_constant_matrix(H, M):
    """docstring for compute_force_constant_matrix"""
    return dot( dot(M, H), M )

def compute_normal_modes(F, discard_trans_and_rot=True):
    eig_vals, eig_vecs = eig(F)
    eig_vals = eig_vals.real
    inds = argsort(eig_vals)

    if discard_trans_and_rot:
        # option to remove translational and rotational modes
        # there are six such modes for for a system with 3N coords
        inds = inds[6:]

    eig_vals = eig_vals[inds]
    eig_vecs = eig_vecs[:,inds]
    modes = NormalModes(eig_vals, eig_vecs)
    return modes

def generate_mode_trajectory(normal_modes, initial_coords, mode_number,
                             peak_scale_factor=1.0, num_steps_per_peak=10):
    zero_to_max = linspace(0.0, peak_scale_factor, num_steps_per_peak)
    max_to_zero = linspace(peak_scale_factor, 0.0, num_steps_per_peak)[1:]
    zero_to_min = linspace(0.0, -peak_scale_factor, num_steps_per_peak)[1:]
    min_to_zero = linspace(-peak_scale_factor, 0.0, num_steps_per_peak)[1:]
    scale_factors = concatenate( [zero_to_max, max_to_zero, zero_to_min,
                                  min_to_zero] )

    trajectory = Trajectory()
    freq, amplitude = normal_modes.get_mode(mode_number)
    assert amplitude.shape[1] == 3, "Cartesian mode amplitude should have" \
                                    "three columns (xyz)."

    for this_scale_factor in scale_factors:
        new_coords = initial_coords * (1 + this_scale_factor * amplitude)
        trajectory.add_frame(new_coords)

    return trajectory


class NormalModes(object):
    """docstring for NormalModes"""
    def __init__(self, eig_vals, eig_vecs, reshape_eig_vecs=True):
        self.frequencies = list(eig_vals)
        self.amplitudes = []
        for i in xrange(eig_vecs.shape[1]):
            this_eig_vec = eig_vecs[:,i]
            if reshape_eig_vecs:
                this_mode_amplitude = this_eig_vec.reshape((len(this_eig_vec)/3,3))
            else:
                this_mode_amplitude = this_eig_vec
            self.amplitudes.append(this_mode_amplitude)

    def __str__(self):
        return "\n%s\n%s" % (self.frequencies, self.amplitudes)

    def get_mode(self, mode_num):
        return self.frequencies[mode_num], self.amplitudes[mode_num]

    def get_frequencies(self):
        return self.frequencies

    def get_amplitudes(self):
        return self.amplitudes

    def freq_to_str(self):
        freq_str = ""
        for i, frequency in enumerate(self.frequencies):
            freq_str += "%d\t%.1f\n" % (i, frequency)
        return freq_str
