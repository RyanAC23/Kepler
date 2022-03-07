"""Script for playing with animated planetary orbits."""

from matplotlib import pyplot as plt
from matplotlib import animation

from ivp_evolvers import Orbits

from planets import *

from planets import Planet

# ----- Class Call --------------------------


class EllipticalOrbit2_10(Planet):
    r"""Object containing parameters for a generic elliptical system.

    Works well with \beta=2.10. Somewhat stable, rapidly precessing.
    Compare with the much lower precession rate for this same data
    set, but with \beta=2.01.
    """

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 1.0  # AU
        self.eccentricity = 1
        self.period = 1  # years

        self.x0 = self.semimajor + 0j
        self.v0 = 1j * 2.5


class EllipticalOrbit3_0(Planet):
    r"""Object containing parameters for a generic elliptical system.

    Works well with \beta=3.0. Elliptical orbits do not close here.
    A perfectly circular orbit will close, but eventually even numerical
    errors will destroy the stability.
    """

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 1.0  # AU
        self.eccentricity = 1
        self.period = 1  # years

        self.x0 = self.semimajor + 0j
        self.v0 = 1j * 5


class EllipticalOrbitMinus4_01(Planet):
    r"""Object containing parameters for a generic elliptical system.

    This one looks nice for these initial conditions and \beta=-4.01.
    """

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 1.0  # AU
        self.eccentricity = 1
        self.period = 1  # years

        self.x0 = self.semimajor + 0j
        self.v0 = 1j * 2.5


MyPlanet = Mercury
c = Orbits(
    planet=MyPlanet, dt=0.001, Tmax=10 * getattr(MyPlanet(), "period"), exponent=2.0
)

# --- static plot ----------------------------

static_plot = False
if static_plot:
    # plt.figure(figsize=(7, 7))
    plt.plot(c.data.real, c.data.imag, "--")
    plt.scatter(0, 0, color="k")
    plt.xlabel("distance (x)")
    plt.ylabel("distance (y)")
    plt.title(rf"Central Potential for $\beta$={c.exponent:.2f}")

    plt.show()

# --- animated plot --------------------------

Lmax = c.x0.real * 2

fig = plt.figure(figsize=(7, 7))
ax = plt.axes(
    xlim=(-Lmax, Lmax), ylim=(-Lmax, Lmax)
)  # can't change this while blitting
(line,) = ax.plot([], [], "k", lw=2)
(points,) = ax.plot([], [], "mo", ms=6)
label = ax.text(
    -Lmax * 0.95, Lmax * 0.80, "", ha="left", va="center", fontsize=16, color="Black"
)


def init():
    """Init data for plot functions."""
    line.set_data([], [])
    points.set_data([], [])
    ax.set_xlabel("distance (x)")
    ax.set_ylabel("distance (y)")
    ax.set_title(rf"Central Potential for $\beta$={c.exponent:.2f}")
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
        f"time={c.ts[i]:.2f}\nenergy={c.Es[i]:.2f}\n"
        + f"KE={c.KEs[i]:.2f}\nPE={(c.Es[i]-c.KEs[i]):.2f}"
    )
    return (line, points, label)


anim = animation.FuncAnimation(
    fig, animate, init_func=init, frames=len(c.ts), interval=1, blit=True
)

plt.show()
