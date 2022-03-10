# Keplerian Orbits

# [Basic Orbits GitHub Website](https://ryanac23.github.io/Kepler/Docs/_build/html/index.html)

Simple implementation of `scipy.solve_ivp` to integrate classical nonrelativistic orbits of solar system planets around the sun. A discussion of the equations of motion is included. At the end, I've worked through an extrapolation problem for computing the perihelion shift of Mercury, which works to a low accuracy.

Examples of the code and how it works are contained in the documentation. The coordinates of the problem are represented by complex numbers, where `x` and `y` are typically represented by `z.real` and `z.imag`, respectively.

### ----- Useful Links -----------------------------------------------------

* [RK4 error in solve_ivp](https://stackoverflow.com/questions/53645649/cannot-get-rk4-to-solve-for-position-of-orbiting-body-in-python#comment94157867_53646267)
* [Set fixed integration step size with solve_ivp](https://stackoverflow.com/questions/54494770/how-to-set-fixed-step-size-with-scipy-integrate)
