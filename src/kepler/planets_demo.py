"""Script for playing with animated planetary orbits."""

from matplotlib import pyplot as plt
from matplotlib import animation

import mmf_setup

mmf_setup.set_path()

from evolvers.ivp_evolvers import Orbits
from kepler.planets import *
from kepler.planets import Planet

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
# ----- Class Call --------------------------
# -------------------------------------------

MyPlanet = Mercury
c = Orbits(
    planet=MyPlanet, dt=0.001, Tmax=10 * getattr(MyPlanet(), "period"), exponent=2.00
)
print(abs(c.data).max(), abs(c.data).min())

# -------------------------------------------
# --- static plot ---------------------------
# -------------------------------------------

static_plot = True
if static_plot:
    # plt.figure(figsize=(7, 7))
    plt.plot(c.data.real, c.data.imag, "--")
    plt.scatter(0, 0, color="k")
    plt.xlabel("distance (x)")
    plt.ylabel("distance (y)")
    plt.title(rf"Central Potential for $\beta$={c.exponent:.2f}")

    plt.show()

# -------------------------------------------
# --- animated plot -------------------------
# -------------------------------------------

Lmax = c.x0.real * 1.6

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
