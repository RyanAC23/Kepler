"""Work for estimating the semimajor precession of Mercury."""
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
from scipy.optimize import curve_fit

import mmf_setup

mmf_setup.set_path()

from evolvers.ivp_evolvers import Orbits
from kepler.planets import Mercury

eps = np.finfo(float).eps

# ----- Global Parameters -------------------------------------------------------

MyPlanet = Mercury

a = np.linspace(1e-4, 9e-4, 9)
b = np.linspace(1e-3, 4e-3, 4)
alphas = np.concatenate([a, b])

dt = 0.001
Tmax = 20 * getattr(MyPlanet(), "period")

# ----- Relativistic Correction Evolver -----------------------------------------


class RelativisticCorrection(Orbits):
    def __init__(self, alpha=0, **kw):
        self.alpha = alpha
        super().__init__(**kw)

    def dq_dt(self, t, q):
        """Return the rhs of the differential equation to integrate."""
        _x, dx = self.unpack(q)

        # Calculate r^2 with absolute value
        r = abs(_x)
        ddx = -self.k * _x / r ** (self.exponent + 1) * (1 + self.alpha / r ** 2)

        return self.pack(dx, ddx)

    def get_energy(self):
        r"""Calculate the energy for a 1/r^2 orbit. Returns [E, KE] as an array.

        Note that this will need to be overridden or recomputed when
        changing \beta.
        """
        KE = 0.5 * abs(self.vs) ** 2
        PE = -self.k / abs(self.data) * (1 + self.alpha / abs(self.data) ** 2 / 3)
        E = KE + PE

        # return E.ravel()
        return np.array([E, KE])[:, :, -1]


# ----- Calculate the precesssion -----

precession_rates = []


def func(t, a, b):
    return a * t + b


def get_precession_rate(MyPlanet, dt, Tmax, alpha):

    c = RelativisticCorrection(planet=MyPlanet, dt=dt, Tmax=Tmax, alpha=alpha)

    peaks = np.where(
        [
            abs(c.data)[j - 1] > abs(c.data)[j] and abs(c.data)[j] < abs(c.data)[j + 1]
            for j in range(1, len(c.data) - 1)
        ]
    )[0]

    fit_times = c.ts[peaks]
    fit_angles = np.arccos(c.data[peaks].real / abs(c.data[peaks])).ravel()

    popt, pcov = curve_fit(
        func, fit_times, fit_angles, bounds=((-np.inf, 0), (np.inf, eps))
    )

    return popt[0]


for i in range(len(alphas)):
    _current_rate = get_precession_rate(MyPlanet, dt, Tmax, alphas[i]) * 180 / np.pi
    precession_rates.append(_current_rate)
    print(f"alpha {i:2d}: {alphas[i]:.2e} AU^2\trate: {_current_rate:.5f} deg/year")


popt_final, pcov_final = curve_fit(
    func, alphas, precession_rates, bounds=((-np.inf, 0), (np.inf, eps))
)

_alphas = np.array(alphas)

mercury_alpha = 1.1e-8  # AU^2
mercury_precession_rate = (
    func(mercury_alpha, popt_final[0], 0) * 3600 * 100
)  # arcseconds / century
print(f"Mercury's Precession Rate : {mercury_precession_rate:.2f} arcseconds / century")
