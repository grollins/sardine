from scipy.optimize import fmin_bfgs
from .trajectory import Trajectory


class BFGSMinimizer(object):
    """
    Attributes
    ----------
    optimization_fcn : callable f(func, x0, *args)
        A function that optimizes parameters based on scores computed
        from a scoring function.

    Parameters
    ----------
    gtol : float, optional
    epsilon : float, optional
    maxiter : int, optional
    """
    def __init__(self, gtol=1e-6, epsilon=0.01, maxiter=200):
        super(BFGSMinimizer, self).__init__()
        self.gtol = gtol
        self.epsilon = epsilon
        self.maxiter = maxiter
        self.traj = None

    def run_minimization(self, energy_fcn, X, num_atoms,
                         save_trajectory=False, noisy=False):
        """
        Optimize parameters based on a scoring function.

        Parameters
        ----------
        energy_fcn : callable f(x, *args)
            A function that computes the energy, given `x`, an array of parameters.
        X : ndarray
            Initial coordinates to pass to the energy function.
        noisy : bool, optional
            Whether to write minimizer messages to stdout.

        Returns
        -------
        parameter_set : ParameterSet
            Minimized position.
        energy : float
            The energy at the minimized position.
        """
        if save_trajectory:
            callback_fcn = self.make_callback_fcn(num_atoms)
            callback_fcn(X) # save initial coords to trajectory

        results = fmin_bfgs(energy_fcn, x0=X, gtol=self.gtol,
                            epsilon=self.epsilon, maxiter=self.maxiter,
                            full_output=noisy, callback=callback_fcn)
        X_min = results[0]
        energy = results[1]
        return X_min, energy

    def make_callback_fcn(self, num_atoms):
        self.traj = Trajectory()
        def callback_fcn(X_vec):
            X = X_vec.reshape((num_atoms, 3))
            self.traj.add_frame(X)
        return callback_fcn

    def get_trajectory(self):
        return self.traj
