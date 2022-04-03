"""Module containing plots and estimates of Halley's Comet's orbit and velocity."""

from matplotlib import pyplot as plt
from matplotlib import animation

import mmf_setup

mmf_setup.set_path()

from evolvers.ivp_evolvers import Orbits
from kepler.planets import HalleysComet

# -------------------------------------------
# ----- Main --------------------------------
# -------------------------------------------

if __name__ == "__main__":
    main()


def main():
    anim = animation.FuncAnimation(
        fig, animate, init_func=init, frames=len(c.ts), interval=1, blit=True
    )

    plt.show()


# -------------------------------------------
# ----- Body --------------------------------
# -------------------------------------------


MyPlanet = HalleysComet
c = Orbits(planet=MyPlanet, dt=0.1, Tmax=3 * getattr(MyPlanet(), "period"))

fig = plt.figure(figsize=(15, 6))
ax = plt.axes(xlim=(-70, 1), ylim=(-10, 10))  # can't change this while blitting
(line,) = ax.plot([], [], "k", lw=2)
(points,) = ax.plot([], [], "mo", ms=6)
label = ax.text(-70 * 0.95, 0, "", ha="left", va="center", fontsize=16, color="Black")
title = ax.set_title("Halley's Comet, orbit=[]")


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
