from collections import defaultdict
import csv
import numpy as np
import pytest

import mmf_setup

mmf_setup.set_path()
from kepler.planets import *
from kepler.planets import get_planet_dict

# For a proper testfile, this should have measurements of all of the values kepler.planets calculates.

planetary_data = get_planet_dict()

planet_names = [
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "Pluto",
]


@pytest.mark.parametrize("planet_name", planet_names)
def test_planet_params(planet_name):
    """Assert checks that the given planet has good data hardcoded or calculated.

    Use of `globals() is frowned upon. Look for a cleaner solution.`"""
    print("Testing %s..." % planet_name)
    P = globals()[planet_name]()

    # assert(np.allclose(P.mass, planetary_data[planet_name]["mass"], atol=1e-1))
    assert np.allclose(P.semimajor, planetary_data[planet_name]["a"], rtol=1e-2)
    assert np.allclose(P.eccentricity, planetary_data[planet_name]["e"], rtol=5e-2)
    assert np.allclose(P.period, planetary_data[planet_name]["p"], rtol=1e-2)
    # assert(np.allclose(P.inclination, planetary_data[planet_name]["I"], atol=1e-1))
    # assert np.allclose(
    #    P.eq_inc, planetary_data[planet_name]["equatorial inclination to orbit"], atol=1e-1
    # )


def test_all_planets():
    for planet in planet_names:

        test_planet_params(planet)


# semimajor_axis[Planet.__name__], getattr(Planet(), "semimajor")
