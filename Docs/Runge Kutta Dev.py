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
#     display_name: Python [conda env:python3]
#     language: python
#     name: conda-env-python3-py
# ---

# For the last time, this notebook will be a dev station to create two Runge Kutta evolver tools.
#
# Make an evolver that does anything first. Euler, who cares. Worry about how to improve later.

# ## Prelim - Euler Method

# $$
# \begin{gather} \
# f(x) = \frac{dx}{dt} \approx \frac{\Delta x}{\Delta t} = \frac{x[i] - x[i-1]}{h} \\ \\
# x[i] = f(x) + x[i-1]h \\
# y(n+1) = y_n + h f(x_n, y_n)
# \end{gather}
# $$

# %pylab inline --no-import-all
import numpy as np


# \begin{equation*}
# \frac{dx}{dt} = f(x, t) = a t,
# \rightarrow x = \frac{at^2}{2}
# \end{equation*}


def euler(func, tf, y0, N):
    # h = dt
    ts = [0.0]
    h = tf / N
    ys = [y0]

    """do all steps now; change to do one step
    with an external looper (evolver) later"""
    for i in range(N):
        y = ys[-1]
        t = ts[-1]
        _y = y + h * func(y0, t)
        ys += [_y]
        ts += [t + h]
        del _y, y, t

    return ys, ts


# +
"""need to define a time step, initial
state and function. Here we will use
dx/dt = v0 + at"""


def func(v0, t):
    a = 9.8  # m/s**2
    return v0 + a * t


N = 10000
v0 = 0.0  # m/s
t_final = 10.0  # s? Final velocity should be ~10 m/s.
# -

xs, ts = euler(func, t_final, v0, N)


# +
def x_an(t):
    return 0.0 + (0.0 * t) + (9.8 * t ** 2 / 2.0)


t_test = np.linspace(0.0, 10.0, N + 1)
x_test = x_an(t_test)

assert np.allclose(ts, t_test)
ts[-1], t_test[-1]
len(ts), len(t_test)
# -

plt.plot(t_test, 0.99 * x_test)
plt.plot(ts, xs)


# This works, though the namespace is very cluttered. Next we'll try to implement a class structure.
#
# Something like
#
# > evolver(euler, RK2, RK4, etc)
#
# > ---- Initialize (physical info common to each evolver)
#
#
# > ---- state (what you actually type out each time)
#
# You should only have to fiddle with state to do a problem, and the other two should be inherited. It's also a good idea to set a default dt=0.01 instead of specifying N, so you can eventually make this a general evolver and not a static one.
#
# Info you need to do a problem:
#
# **dt**
#
# **t0**
#
# **tf**
#
# **xyz** (position for dynamics, all points on a grid for a field)

# +
class Evolve(object):
    def __init__(dt=0.01):
        pass

    def euler(func, tf, y0, N):
        pass

    # h = dt
    ts = [0.0]
    h = tf / N
    ys = [y0]

    """do all steps now; change to do one step
    with an external looper (evolver) later"""
    for i in range(N):
        y = ys[-1]
        t = ts[-1]
        _y = y + h * func(y0, t)
        ys += [_y]
        ts += [t + h]
        del _y, y, t

    return ys, ts


class State(object):
    pass


# -


# ## Runge Kutta 2nd Order


# ## Runge Kutta 4th Order


# ### pytimeode example
# \begin{equation*} \frac{\mathrm{d}y(t)}{\mathrm{d}t} = f\bigl(y(t), t\bigr) = -y^2, \qquad y(0) = y_0 = \begin{pmatrix} 1\\ 2 \end{pmatrix} \end{equation*}
#
