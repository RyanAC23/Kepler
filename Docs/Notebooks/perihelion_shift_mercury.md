---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.7
kernelspec:
  display_name: Python 3 (kepler)
  language: python
  name: kepler
---

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

import mmf_setup; mmf_setup.set_path()

from kepler.planets import Mercury
from evolvers.ivp_evolvers import Orbits

eps = np.finfo(float).eps
```

# Mercury

Here I'll take a look at Mercury, which has a noticeable eccentricity, though less extreme than Halley's Comet. In the real world, this orbit precesses, but in simulation it will not. Why the difference?

```{code-cell} ipython3
MyPlanet = Mercury

c = Orbits(planet=MyPlanet, dt=0.0001, Tmax=1.5 * getattr(MyPlanet(), "period"))

from IPython.display import clear_output

fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlabel("distance (AU)")
ax.set_ylabel("distance (AU)")
ax.set_title("Mercury's Orbit")
ax.set(xlim=(-0.6, 0.5), ylim=(-0.5, 0.5), aspect=1)
(line,) = ax.plot([], [], "k", lw=2)
(points,) = ax.plot([], [], "mo", ms=6)

frames = 100
for i in range(len(c.data[::frames])):
    j = int(i * frames)
    x = c.data.real[:j, 0]
    y = c.data.imag[:j, 0]
    line.set_data(x, y)
    points.set_data([0, c.data.real[j, 0]], [0, c.data.imag[j, 0]])

    display(fig)
    clear_output(wait=True)
```

## Deviation from $1/r^2$
In our universe, there are more gravitating bodies than just the sun and Mercury. The contributions to the gravitational force from these other bodies act as a perturbation on Mercury's orbit and cause it to precess. One way of roughly modelling this is by slightly deviating from the $1/r^2$ force.

```{code-cell} ipython3
MyPlanet = Mercury

c = Orbits(
    planet=MyPlanet, dt=0.0001, Tmax=2 * getattr(MyPlanet(), "period"), exponent=2.2
)

from IPython.display import clear_output

fig, ax = plt.subplots(figsize=(7, 7))
ax.set_xlabel("distance (AU)")
ax.set_ylabel("distance (AU)")
ax.set_title(fr"Mercury's Orbit, $\beta$={c.exponent:.2f}")
ax.set(xlim=(-0.4, 0.4), ylim=(-0.4, 0.4), aspect=1)
(line,) = ax.plot([], [], "k", lw=2)
(points,) = ax.plot([], [], "mo", ms=6)

frames = 100
for i in range(len(c.data[::frames])):
    j = int(i * frames)
    x = c.data.real[:j, 0]
    y = c.data.imag[:j, 0]
    line.set_data(x, y)
    points.set_data([0, c.data.real[j, 0]], [0, c.data.imag[j, 0]])

    display(fig)
    clear_output(wait=True)
```

# Measuring the perihelion shift

$$
\frac{\partial^2 r}{\partial t^2} = \frac{k}{r^2}\bigg( 1 + \frac{\alpha}{r^2}\bigg)
$$

A quick note - the aphelion could have been used as easily as the perihelion, right? Not quite. Since the precession (the way the problem has been defined) proceeds counter-clockwise, the angle of perihelion should be increasing in time. Since $\arccos$ has a range of $[0,\pi]$, if the angle is more than $\pi$, the arguments of $\arccos$ will give the angle as reflected across the x-axis. So instead of $\pi + \delta$, we would measure $\pi - \delta$ rad. This problem goes away if we use the perihelion, which starts at angle $0$ rad.

```{code-cell} ipython3
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
```

## Precession rate curve fit test

```{code-cell} ipython3
alphas = [6e-3, 5e-3, 3e-3, 1e-3, 5e-4]
alphas = [alphas[-1]]
precession_rates = []


def func(t, a, b):
    return a * t + b


for i in range(len(alphas)):

    MyPlanet = Mercury
    c = RelativisticCorrection(
        planet=MyPlanet,
        dt=0.0001,
        Tmax=10 * getattr(MyPlanet(), "period"),
        alpha=alphas[i],
    )

    peaks = np.where(
        [
            abs(c.data)[i - 1] > abs(c.data)[i] and abs(c.data)[i] < abs(c.data)[i + 1]
            for i in range(1, len(c.data) - 1)
        ]
    )[0]

    fit_times = c.ts[peaks]
    fit_angles = np.arccos(c.data[peaks].real / abs(c.data[peaks])).ravel()

    popt, pcov = curve_fit(
        func, fit_times, fit_angles, bounds=((-np.inf, 0), (np.inf, 0.00000001))
    )
    precession_rates.append(popt[0])
```

```{code-cell} ipython3
do_static_plots = True

fit_locations = abs(c.data[peaks]).ravel()
fit_angles = np.arccos(c.data[peaks].real / abs(c.data[peaks])).ravel()

if do_static_plots:
    plt.figure(figsize=(12, 12))
    plt.subplot(221)
    plt.plot(fit_times, fit_angles, ".")
    plt.plot(fit_times, func(fit_times, *popt), "--")

    plt.title(rf"angle of perihelion vs time, $\dot\theta$={popt[0]:.3f} rad / yr")
    plt.xlabel("time (years)")
    plt.ylabel("Angle of max radius")

    plt.subplot(222)
    plt.plot(c.data.real, c.data.imag, "--")
    plt.scatter(0, 0, color="k")
    plt.xlabel("distance (x)")
    plt.ylabel("distance (y)")
    plt.suptitle(rf"Precession of Mercury, $\alpha$ = {c.alpha:.1e}")
    plt.title("Orbital shift")

    plt.subplot(223)
    plt.plot(c.ts, abs(c.data))
    plt.plot(c.ts[peaks], fit_locations, ".")
```

## Extrapolate precession from changing $\alpha$
The above appears to have given a perihelion shift. We now repeat this for many values of $\alpha$, and curve-fit the resulting $\omega$ vs $\alpha$ plot to extrapolate the real value.

Watch out for values of alpha that are too large.

```{code-cell} ipython3
# alphas = [1e-4, 5e-5, 1e-4, 5e-4, 1e-3, 2e-3, 3e-3, 4e-3]

a = np.linspace(1e-4, 9e-4, 9)
b = np.linspace(1e-3, 4e-3, 4)
alphas = np.concatenate([a, b])

precession_rates = []


def func(t, a, b):
    return a * t + b


for i in range(len(alphas)):

    MyPlanet = Mercury
    c = RelativisticCorrection(
        planet=MyPlanet,
        dt=0.001,
        Tmax=20 * getattr(MyPlanet(), "period"),
        alpha=alphas[i],
    )

    peaks = np.where(
        [
            abs(c.data)[i - 1] > abs(c.data)[i] and abs(c.data)[i] < abs(c.data)[i + 1]
            for i in range(1, len(c.data) - 1)
        ]
    )[0]

    fit_times = c.ts[peaks]
    fit_locations = abs(c.data[peaks]).ravel()
    fit_angles = np.arccos(c.data[peaks].real / abs(c.data[peaks])).ravel()

    popt, pcov = curve_fit(
        func, fit_times, fit_angles, bounds=((-np.inf, 0), (np.inf, eps))
    )

    precession_rates.append(popt[0])

    print(popt)
```

```{code-cell} ipython3
precession_rates_degrees = np.array(precession_rates) * 180 / np.pi


def func(t, a, b):
    return a * t + b


popt_final, pcov_final = curve_fit(
    func, alphas, precession_rates_degrees, bounds=((-np.inf, 0), (np.inf, eps))
)

_alphas = np.array(alphas)

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(alphas, precession_rates_degrees, ".")
ax.plot(alphas, func(_alphas, *popt_final), "--")
ax.set(xlabel=r"$\alpha$", ylabel=r"$\omega$ (rad / year)")
ax.set(title=fr"Precession rate vs $\alpha$: slope={popt_final[0]:.3f}")
plt.tight_layout()


mercury_alpha = 1.1e-8  # AU^2
mercury_precession_rate = (
    func(mercury_alpha, popt_final[0], 0) * 3600 * 100
)  # arcseconds / century
print(f"Mercury's Precession Rate : {mercury_precession_rate:.2f} arcseconds / century")
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```
