"""Simple evolvers for quick implementation of IVP problems."""

import numpy as np
from scipy.integrate import solve_ivp
import warnings


# ----- Simple Gravitational Orbits ---------------------------------------


class Orbits:
    """Class solving for the trajectories of central potential orbits."""

    def __init__(self, planet=None, dt=0.1, Tmax=10, k=None, exponent=2):
        """Calculate 2D trajectories in a central potential problem.

        parameters
        ---------------------
        k: float
            Gravitational coupling constant
        exponent: float
            Factor n in 1/r^n in force equation
        """
        if k is None:
            self.k = 4 * np.pi ** 2  # AU^3 / yr^2
        else:
            self.k = k

        self.exponent = exponent
        self.dt = dt
        self.Tmax = Tmax
        self.ts = np.arange(0, self.Tmax, self.dt)

        if planet is None:
            warnings.warn("No planet data supplied. Using default values for Earth.")
            self.x0 = 1.0 + 0j
            self.v0 = 0.0 + 2 * np.pi * 1j
        else:
            self.x0 = getattr(planet(), "x0")
            self.v0 = getattr(planet(), "v0")

        xs, vs = self.integrate()
        self.data = xs
        self.vs = vs
        self.Es, self.KEs = self.get_energy()

    def initial_state(self):
        """Return initial conditions for IVP."""
        q0 = self.x0
        dq0 = self.v0

        return [q0, dq0]

    def pack(self, x, dx):
        """Ravel two arrays into one for passing into an integrator."""
        return np.ravel([x, dx])

    def unpack(self, q):
        """Turn one packed array into two unpacked ones."""
        new_shape = (2, q.shape[0] // 2) + q.shape[1:]
        return q.reshape(new_shape)

    def dq_dt(self, t, q):
        """Return the rhs of the differential equation to integrate."""
        _x, dx = self.unpack(q)

        # Calculate r^2 with components
        x, y = _x.real, _x.imag
        r2 = (x ** 2 + y ** 2) ** ((self.exponent + 1) / 2)
        _ddx = -self.k * (x + 1j * y) / r2

        # Calculate r^2 with absolute value
        r = abs(_x)
        ddx = -self.k * _x / r ** (self.exponent + 1)

        assert np.allclose(_ddx, ddx)

        return self.pack(dx, ddx)

    def integrate(self):
        """Solve for the orbital trajectories.

        Note that solve_ivp has issues with these orbital problems when
        using its adaptive step size integration methods. To ensure
        convergence, rtol and atol have been set here to be low enough that
        issues with the adaptive step size don't appear (the energy of the
        orbiting body does not change in systems where energy should be
        conserved), but also as high as possible to keep adaptive measures
        possible. If the energy changes for a 1/r potential, consider lowering
        both atol and rtol further.
        """
        q0 = self.initial_state()

        res = solve_ivp(
            self.dq_dt,
            [0, self.Tmax],
            q0,
            t_eval=self.ts,
            method="BDF",
            dense_output=True,
            atol=0.000001,  # to force non-adaptive step sizes
            rtol=0.000001,
            max_step=self.dt,
        )
        ts, (zs, dzs) = res.t, self.unpack(res.y)

        self.ts = ts

        zs = zs.T
        dzs = dzs.T
        return zs, dzs

    def get_energy(self):
        r"""Calculate the energy for a 1/r^2 orbit. Returns [E, KE] as an array.

        Note that this will need to be overridden or recomputed when
        changing \beta.
        """
        KE = 0.5 * abs(self.vs) ** 2
        PE = -self.k / abs(self.data)
        # PE = -self.k / (self.exponent - 1) / abs(self.data) ** (self.exponent - 1)
        # figure out which cases this won't apply for and flag them
        E = KE + PE

        # return E.ravel()
        return np.array([E, KE])[:, :, -1]
