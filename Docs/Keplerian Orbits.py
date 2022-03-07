# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python [conda env:work]
#     language: python
#     name: conda-env-work-py
# ---

# # Gravitational Orbits

# +
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

from matplotlib import pyplot as plt
from matplotlib import animation

# -

# ## Equations of motion
#
# For an inverse square attractive potential, the energy takes the form
#
# $$
# E = \frac{1}{2} m v^2 - \frac{km}{r} \qquad k \equiv GM
# $$
# In the reduced mass picture,
# $$
# r = |\vec{r_1} - \vec{r_2}| \qquad \frac{m_1 m_2}{m_1 + m_2} = \mu = m
# $$
# Since for most of the work in this notebook the larger planet's mass (i.e., the sun) is hidden in the coupling constant $k$, I won't make an effort to use $\mu$ instead of $m$. This would be important if I was looking at three body interactions or two body orbits where the masses are comparable and the center of mass is outside of either body, but that won't be the case here. Be aware that it could cause trouble in the future and that this code will not be completely generic.
#
# ### Lagrangian
#
# Though trivial up to a certain point, it's worth refreshing myself on the construction of the problem from scratch.
#
# $$
# L = \frac{1}{2}m|\dot{\vec{r}}|^2 - U(r) \qquad U(r) = - \frac{km}{r} \\
# \frac{1}{2}m|\dot{\vec{r}}|^2 = \frac{1}{2}m |\vec{v_x} + \vec{v_y}|^2 \\
# x = r \cos \phi \qquad y = r \sin \phi \\
# \dot x =  \dot r \cos \phi - \dot \phi r \sin \phi \qquad \dot y =  \dot r \sin \phi + \dot \phi r \cos \phi \\
# \rightarrow L = \frac{1}{2}m(\dot r \cos \phi \hat x - \dot \phi r \sin \phi \hat x + \dot r \sin \phi \hat y + \dot \phi r \cos \phi \hat y)^2 - U(r)\\
# = \frac{1}{2} (\dot \phi ^2 r ^2 + \dot r ^2 ) - U(r)
# $$
#
# Because the Lagrangian depends only on the radius, the angular momentum will be conserved. Putting these through the Euler-Lagrange equations,
#
# $$
# \frac{\partial}{\partial t} \frac{\partial L}{\partial \dot r} = \frac{\partial L}{\partial r} \qquad \frac{\partial}{\partial t} \frac{\partial L}{\partial \dot \phi} = \frac{\partial L}{\partial \phi} = 0 \\
# m \ddot r = m r \dot \phi ^2 - \frac{\partial U(r)}{\partial r} \qquad \partial_t (\dot \phi m r^2) = 0 \\
# \rightarrow \dot \phi = \frac{\mathcal{l}}{mr^2}
# $$

# The $r$ equation gives us the force we could have derived from Newton's laws, and the $\phi$ equation shows us that the angular momentum is a symmetry of the system. Because of this, our orbits will be constrained to a 2D plane.
#
# ## Effective potential
#
# Putting these together yields the force equation
#
# $$
# F_r = \underbrace{\frac{l^2}{m r^3} - \frac{\partial U(r)}{\partial r}}_{U_\text{eff}}
# $$
#
# Calling the angular momentum piece an effective potential $U_\text{eff}$ is common, as is seeing the radial potental, the angular momentum term, and $U_\text{eff}$ all plotted together for demonstration. I'll do that below for an inverse square potential (I haven't used a specific form of the potential in the above, aside from constraining it to depend only on r!) with $k=1$.

# +
k = 24 / 100
m = 1 / 100
l = 1.3 / 100

rs = np.linspace(0.04, 0.3, 100)


def inverse_square(r):
    return -k / r


def angular_term(r):
    return l ** 2 / 2 / m / r ** 2


def U_effective(r):
    return angular_term(r) + inverse_square(r)


plt.figure(figsize=(5, 5))
plt.xlim(0.025, 0.2)
plt.ylim(-4, 4)
plt.plot(rs, inverse_square(rs), label=r"$-k/r$")
plt.plot(rs, angular_term(rs), label=r"$l^2/2mr^2$")
plt.plot(rs, U_effective(rs), label=r"$U_{eff}$")
plt.legend()
plt.tight_layout()


# + [markdown] slideshow={"slide_type": "slide"}
# It turns out it's a bit tough to find good parameters to demonstrate the potential curve. The green curve is the effective potential; bound orbits are when the energy of the orbiting body is in the concave region, below the region where it would scatter (unbound orbit).
#
# # 1D Central Potential Scattering
#
# The total energy is
#
# $$
# E = \frac{1}{2} m \dot r^2  + U_{\text{eff}} \\
# \frac{\partial r}{\partial t} = \sqrt{\frac{2}{m} \big( E - U_{\text{eff}} \big)} \\
# \int_0^t dt = \int_0^R \frac{dr}{\sqrt{\frac{2}{m} \big( E - U_{\text{eff}} \big)}}
# $$
#
# This is where the work stops being trivial, because these integrals are hard. The time as a function of the energy and radius is also a somewhat bogus thing to solve for. Inverting this equation does not look like a fun time. The trick (which makes things easier but not necessarily EASY) is to change coordinates as follows:
#
# $$
# \frac{\partial r}{\partial t} = \frac{\partial r}{\partial \phi}\frac{\partial \phi}{\partial t} =
# \frac{\partial r}{\partial \phi} \frac{\mathcal l}{mr^2}
# $$
#
# Then
#
# $$
# E = \frac{1}{2}\frac{\mathcal l^2}{mr^4} \bigg(\frac{\partial r}{\partial \phi}\bigg)^2 + U_{\text{eff}}
# $$
#
# $$
# \boxed{\int_{\phi_0}^\phi = \Delta \phi = \int_0^R \frac{dr}{r^2 \sqrt{\frac{2m}{\mathcal l^2} \big( E - U_{\text{eff}} \big)}}}
# $$

# + [markdown] slideshow={"slide_type": "slide"}
# For the inverse square problem,
#
# $$
# \Delta \phi = \int_0^R \frac{dr}{r^2 \sqrt{\frac{2m}{\mathcal l^2} \big( E + \frac{k}{r} - \frac{\mathcal l^2}{2mr^2} \big)}} \\
# = \int_0^R \frac{dr}{r^2 \sqrt{\big(\frac{2mE}{l^2} + \frac{2mk}{l^2r} - \frac{1}{r^2} \big)}}\\
# = - \int \frac{du}{\sqrt{\big(\frac{2mE}{l^2} + \frac{2mk u}{l^2} - u \big)}}
# $$
#
# The last line comes from making the substitution $u=1/r$, $du = -1/r^2 dr$.
#
# ## Elliptic integral gratis
# This is gross. I don't like it. I'm not good at these things. Here's the integral:
# $$
# \int \frac{dx}{\sqrt{\alpha + \beta x + \gamma x^2}} = \frac{1}{\sqrt{-\gamma}} \arccos \bigg( - \frac{\beta + 2 \gamma x}{\sqrt{\beta^2 - 4 \alpha \gamma}}\bigg)
# $$
#
# It looks like some completing the square action is going on, given the $\arccos$ argument in the denominator.
#
# Making the substitutions and using the integral from the table....
#
# $$
# alpha = \frac{2mE}{l^2} \qquad \beta = \frac{2mk}{l^2} \qquad \gamma = -1 \\
# \Delta \phi = \arccos\bigg(\frac{\frac{l^2 u}{2mk} - 1}{\sqrt{1 + \frac{2El^2}{mk^2}}}\bigg)\\
# $$
#
# Solving for $u$ and returning to $u = 1/r$,
#
# $$
# \boxed{\frac{1}{r} = \bigg(\frac{2mk}{l^2}\bigg)\bigg(1 + \sqrt{1 + \frac{2El^2}{mk^2}}\cos\Delta\phi \bigg)}
# $$
# Smart people realized that this is the equation for conic sections,
#
# $$
# \frac{1}{r} = C(1 + \epsilon \cos \Delta \phi ) \\
# C = \bigg(\frac{2mk}{l^2}\bigg) \qquad \epsilon = \sqrt{1 + \frac{2El^2}{mk^2}}
# $$
#
#
# * Physical insight to derive these damned equations #1:
#
# **At the farthest and nearest points, the velocity of the orbiting body is entirely tangential. At these points, the radial velocity is zero.**
#     $$
#     E = \frac{-k}{r} + \frac{l^2}{2mr^2} \\
#     E + \frac{k}{r} - \frac{l^2}{2mr^2} = 0 \\
#     \rightarrow r = \frac{-k \pm \sqrt{k^2 + \frac{2El^2}{m}}}{2E}
#     $$
#     Using the definition of $\epsilon$ this is
#     $$
#     r = (1 \pm \epsilon) \qquad a \equiv \frac{k}{2E}\\
#     r_a = a(1+\epsilon) \qquad r_p = a(1-\epsilon)
#     $$
#     $r_a$ is the radius of the orbit at apoapsis, and $r_b$ is the radius of the orbit at periapsis (when discussing orbits around our sun, these terms take on the more specific *perihelion* and *aphelion*). $a$ is the semimajor axis. Note that $r_a + r_p = 2a$.
#
# ---
#
# There's one more thing that I need to get reasonable simulations for the objects in our solar system. Somehow I need to estimate the initial velocity. The easiest way I could figure is when the velocity is entirely in the tangential direction, at one of the apsides. In these cases we don't need to guess, we can calculate.
#
# [a picture you've drawn out of the system would be nice here.]
#
# If $b$ is the semiminor axis, we can relate the angular momentum (conserved!) at any point with that at b.
# $$
# l_1 = l_2 \qquad m v_1 r_1 = m v_2 r_2 = m b v_2
# $$
# The mass drops out, and the relationship is an easy one. Combine this with conservation of energy between the two points:
#
# $$
# \frac{1}{2} v_1 ^2 - \frac{k}{r_1} = \frac{1}{2} v_2 ^2 - \frac{k}{b}
# $$
#
# Two things of use:
# $$
# r_{\text{apsis}} = a(1\pm \epsilon) \qquad b = \sqrt{r_a r_b} \rightarrow b = a\sqrt{1-\epsilon^2}
# $$
#
# By combining the above equations and doing a few pages of dumb messy algebra (I did it, it works!) the velocities are
#
# $$
# v_a = \sqrt{
#             \frac{k}{a} \frac{(1-\epsilon)}{(1+\epsilon)}
#             }
#       \qquad
# v_p = \sqrt{
#             \frac{k}{a} \frac{(1+\epsilon)}{(1-\epsilon)}
#             }
# $$

# + [markdown] slideshow={"slide_type": "slide"}
# # Gravitational Orbits Code
#
# Using the above equations, everything necessary to start simulating orbits is done. There are other things to discuss, but this has been enough equations for now. Let's break it up with the main code for the orbital routine.
# -


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
            # max_step=self.dt,
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
        PE = -self.k / (self.exponent - 1) / abs(self.data) ** (self.exponent - 1)
        E = KE + PE

        # return E.ravel()
        return np.array([E, KE])[:, :, -1]


# This gives us most of the necessary tools to simulate orbits for a general exponent $\beta$ in a $k / r^{\beta}$ force law. To be clear, in our universe (as far as we know), this exponent is $\beta = 2$. We need planetary data, which is provided below. The attributes specific to the planet of interest get pulled into the `Orbits` class when properly applied.

# ## Define celestial objects

# +
class Planet(object):
    """Example class for Planet objects."""

    k = 4 * np.pi ** 2  # in AU and years

    def __init__(self):
        """Necessary simulation parameters for orbits.

        parameters
        -------------------------
        mass : float
            The mass of the planet.
        semimajor : float
            Length of the semimajor axis of the planet in AU.
        eccentricity : float
            A measure of the deviation from a circular orbit.
        x0 : complex
            The initial position in x + iy.
        v0 : complex
            Initial velocity in x + iy.
        """
        self.mass = None
        self.semimajor = None  # AU
        self.eccentricity = None
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = None
        self.v0 = self.velocity_at_perihelion()

    def velocity_at_perihelion(self):
        """Calculates the velocity at the perihelion."""
        v_p = np.sqrt(
            (self.k / self.semimajor)
            * (1 + self.eccentricity)
            / (1 - self.eccentricity)
        )
        return v_p * 1j


class HalleysComet(Planet):
    """Object containing parameters for Halley's Comet."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 17.955  # AU
        self.eccentricity = 0.96714
        self.period = 76  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.v0 = self.velocity_at_perihelion() * 1j

        self.x0 = self.perihelion + 0j
        #        # self.v0 = 1j * 11.4727
        self.v0 = self.velocity_at_perihelion()


class Mercury(Planet):
    """Object containing parameters for Mercury."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 0.39  # AU
        self.eccentricity = 0.206
        self.period = 87.97 / 365.26  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        # self.v0 = 0 + 2 * np.pi * np.abs(self.x0) / self.period * 1j  # v = gm
        self.v0 = self.velocity_at_perihelion()


class Earth(Planet):
    """Object containing parameters for Earth."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 1.0  # AU
        self.eccentricity = 0.017
        self.period = 1.0  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()


# -

# We don't really need the planet's mass for the problems in this notebook, since the sun is so massive that their total mass is effectively the sun's mass and the center of mass is deep inside that star.
#
# ## units
#
# Let's talk about units.
#
# The orders of magnitude in astronomical work are too spread to be smart for computational work with numbers stored in floats and whatnot. $10^{27}\text{kg}$ doesn't play well with the gravitational constant $G = 6.673 \times 10^{-11}\text{N/kg m}^2$. For problems dealing with our solar system, timescales of orbits can be described between months and hundreds of years, so measuring time in years is safe. We can also measure distance in astronomical units - Pluto is only 39 AU from the sun on average, so we can comfortably go to 100 AU and not push the limits of too many decades in magnitude.
#
# The velocity of a circular orbit is given by Newton's equations:
#
# $$
# \frac{m v^2}{r} = F_r = \frac{GMm}{r^2} \\
# v^2 r = GM
# $$
#
# The Earth's orbit is roughly circular, and it orbits $2\pi \times 1AU$ every year.
#
# $$
# 4 \pi^2 \frac{\text{AU}^2}{1 \text{yr}^2} \times 1 \text{AU} = GM \\
# GM = 4 \pi^2 \frac{\text{AU}^3}{\text{yr}^2} \equiv k
# $$

# ## Earth

# +
MyPlanet = Earth

c = Orbits(planet=MyPlanet, dt=0.0001, Tmax=1.5 * getattr(MyPlanet(), "period"))

from IPython.display import clear_output

fig, ax = plt.subplots()
ax.set_xlabel("distance (AU)")
ax.set_ylabel("distance (AU)")
ax.set_title("Earth's Orbit")
ax.set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), aspect=1)
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
# -

# ## Halley's Comet
# We can check this for Halley's Comet for a more complicated example. The period should be about 76 years.

# +
MyPlanet = HalleysComet

c = Orbits(planet=MyPlanet, dt=0.01, Tmax=1.5 * getattr(MyPlanet(), "period"))

from IPython.display import clear_output

fig, ax = plt.subplots(figsize=(12, 4))
ax.set_xlabel("distance (AU)")
ax.set_ylabel("distance (AU)")
ax.set(xlim=(-40, 5), ylim=(-10, 10), aspect=1)
(line,) = ax.plot([], [], "k", lw=2)
(points,) = ax.plot([], [], "mo", ms=6)

frames = 100
for i in range(len(c.data[::frames])):
    j = int(i * frames)
    x = c.data.real[:j, 0]
    y = c.data.imag[:j, 0]
    ax.set_title(f"Halley's Comet, T={c.ts[j]:.2f} years")
    line.set_data(x, y)
    points.set_data([0, c.data.real[j, 0]], [0, c.data.imag[j, 0]])

    display(fig)
    clear_output(wait=True)
# -

# # Mercury
#
# Here I'll take a look at Mercury, which has a noticeable eccentricity, though less extreme than Halley's Comet. In the real world, this orbit precesses, but in simulation it will not. Why the difference?

# +
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
# -

# ## Deviation from $1/r^2$
# In our universe, there are more gravitating bodies than just the sun and Mercury. The contributions to the gravitational force from these other bodies act as a perturbation on Mercury's orbit and cause it to precess. One way of roughly modelling this is by slightly deviating from the $1/r^2$ force.

# +
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
# -

# # P.S. - Matplotlib AnimateFunc
# The animated plots above look a lot nicer outside of Jupyter notebooks, when using Matplotlib's `AnimateFunc`. Try running one through a terminal!

# +
# %matplotlib notebook

MyPlanet = HalleysComet
c = Orbits(planet=MyPlanet, dt=0.1, Tmax=10 * getattr(MyPlanet(), "period"))

fig = plt.figure(figsize=(8, 4))
ax = plt.axes(xlim=(-70, 1), ylim=(-10, 10))  # can't change this while blitting
(line,) = ax.plot([], [], "k", lw=2)
(points,) = ax.plot([], [], "mo", ms=6)
label = ax.text(-70 * 0.95, 0, "", ha="left", va="center", fontsize=16, color="Black")


def init():
    """Init data for plot functions."""
    line.set_data([], [])
    points.set_data([], [])
    ax.set_xlabel("distance (AU)")
    ax.set_ylabel("distance (AU)")
    ax.set_title("Halley's Comet orbit, T=76 years")
    ax.set(aspect=1)
    return (line, points, label)


def animate(i):
    """Set data for plot functions while looping over index i."""
    x = c.data.real[:i, 0]
    y = c.data.imag[:i, 0]
    line.set_data(x, y)
    line.set_linestyle("--")
    points.set_data([0, c.data.real[i, 0]], [0, c.data.imag[i, 0]])
    label.set_text(
        f"time={c.ts[i]:.2f} years\nenergy={c.Es[i]:.2f}\n"
        + f"KE={c.KEs[i]:.2f}\nPE={(c.Es[i]-c.KEs[i]):.2f}"
    )
    return (line, points, label)


anim = animation.FuncAnimation(
    fig, animate, init_func=init, frames=len(c.ts), interval=1, blit=True
)

# -

# Note that these `Matplotlib.FuncAnimation` plots work much better outside of Jupyter notebooks.
