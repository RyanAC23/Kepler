---
jupyter:
  jupytext:
    formats: ipynb,py:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.3.3
  kernelspec:
    display_name: Python [conda env:work]
    language: python
    name: conda-env-work-py
---

<!-- #region toc=true -->
<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Variations-of-Finite-Difference" data-toc-modified-id="Variations-of-Finite-Difference-1">Variations of Finite Difference</a></span><ul class="toc-item"><li><span><a href="#Euler-Method" data-toc-modified-id="Euler-Method-1.1">Euler Method</a></span></li><li><span><a href="#Runge-Kutta-2nd-Order" data-toc-modified-id="Runge-Kutta-2nd-Order-1.2">Runge-Kutta 2nd Order</a></span></li><li><span><a href="#Runge-Kutta-4th-Order" data-toc-modified-id="Runge-Kutta-4th-Order-1.3">Runge-Kutta 4th Order</a></span><ul class="toc-item"><li><span><a href="#Ex.-Decay" data-toc-modified-id="Ex.-Decay-1.3.1">Ex. Decay</a></span></li><li><span><a href="#Ex.-2nd-order" data-toc-modified-id="Ex.-2nd-order-1.3.2">Ex. 2nd order</a></span></li></ul></li></ul></li></ul></div>
<!-- #endregion -->

# Variations of Finite Difference


## Euler Method


$$y_{n+1} = y_n + h f(x_n, y_n) $$

```python
%pylab inline --no-import-all
```


```python
def euler(f, x0, t0, tf, N):
    h = (tf - t0) / N
    xs = [0]
    ts = [0]
    xs[0] = x0
    ts[0] = t0

    for i in range(N):
        x = xs[-1]
        t = ts[-1]
        x_temp = x + h * f(x, t)
        # How to un-grossify the variables in f above? so that the
        # Func just looks like f(x, t)
        xs += [x_temp]
        ts += [t + h]
        del x_temp

    return ts, xs
```

```python
def f(x, t):
    return t ** 3 + t ** 2 + t + 100
```


## Runge-Kutta 2nd Order


```python
def rk2(f, x0, t0, tf, N):
    h = (tf - t0) / N
    xs = [0]
    ts = [0]
    xs[0] = x = x0
    ts[0] = t = t0

    for i in range(N):
        x = xs[-1]
        t = ts[-1]
        k = [None, None]
        k[0] = h * (f(x, t))
        k[1] = h * (f(x + 0.5 * k[0], t=t + h / 2))
        x_temp = xs[-1] + k[1]
        xs += [x_temp]
        ts += [ts[-1] + h]
        del x_temp, k

    return ts, xs
```


## Runge-Kutta 4th Order


```python
def rk4(f, x0, t0, tf, N):
    h = (tf - t0) / N
    xs = [0]
    ts = [0]
    xs[0] = x = x0
    ts[0] = t = t0

    for i in range(N):
        x = xs[-1]
        t = ts[-1]
        k = [None, None, None, None]
        k[0] = h * (f(x, t))
        k[1] = h * (f(x + 0.5 * k[0], t=t + h / 2.0))
        k[2] = h * (f(x + 0.5 * k[1], t=t + h / 2.0))
        k[3] = h * (f(x + k[2], t=t + h))
        x_temp = xs[-1] + (k[0] + 2 * k[1] + 2 * k[2] + k[3]) / 6.0
        xs += [x_temp]
        ts += [ts[-1] + h]
        del x_temp, k

    return ts, xs
```

```python
plt.figure(figsize=(10, 10))
ts, xs = euler(f, 0, -10, 10, 10)
# have to take much smaller step size with euler to make it match
eu = plt.plot(ts, xs)

ts2, xs2 = rk2(f, 0, -10, 10, 10)
ru2 = plt.plot(ts2, xs2)


ts4, xs4 = rk4(f, 0, -10, 10, 10)
plt.plot(ts4, xs4)
plt.xlabel("ts")
plt.ylabel("xs")

plt.plot(ts, f(xs, ts))
```

### Ex. Decay
$$\frac{dN}{dt} = -\frac{N}{\tau}$$

$$\int_{N_0}^{y}\frac{dN}{N} =\int_{t=0}^{t} -\frac{dt}{\tau}$$

$$ ln\bigg(\frac{N}{N_0}\bigg) = -\frac{t}{\tau}$$

$$ N(t) = N_0 e^{-t / \tau} $$


```python
def dNdt(N, t):
    tau = 2.0
    return -N / tau
```

```python
ts, Ns = euler(dNdt, 100, 0, 10, 10)
ts2, Ns2 = rk2(dNdt, 100, 0, 10, 10)
ts4, Ns4 = rk4(dNdt, 100, 0, 10, 10)
plt.plot(ts, Ns)
plt.plot(ts2, Ns2)
plt.plot(ts4, Ns4)
plt.xlabel("time")
plt.ylabel("number of particles")


# How best to make this (class) inherit vars N0 and tau?
# Would need to map ts to N(t) function, changed in py3
# def N_analytic(t):
#    N0 = 100
#    tau = 2.
#    Ns = []
#    for i in range(len(t)):
#        Ns += N0 * np.exp(-1*t/tau)
#    return Ns
# N_analytics = N_analytic(ts)
# plt.plot(ts, N_analytics)
```

### Ex. 2nd order


```python
def ddt(q, t):
    x, dx = q
    ddx = a = -9.81
    dx = v = a * t
    dq = (dx, ddx)
    return dq
```

```python
%debug
```

```python
ts, dxs = euler(f=ddt, x0=[10, 0], t0=0, tf=10, N=100)
```
