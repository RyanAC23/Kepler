"""Module containing our Solar System's planets and some of their parameters."""
import mmf_setup

mmf_setup.set_path()

import numpy as np

__all__ = [
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "Pluto",
    "HalleysComet",
]

import os


def get_planet_dict():
    planet_data_dict = {}
    abs_path = os.path.dirname(__file__)
    file_path = abs_path + "/planets.csv"

    with open(file_path, "r") as data_file:
        row1 = data_file.readline()
        row1 = row1.strip().split(",")[1:]  # throw away first label
        for row in data_file:
            row = row.strip().split(",")
            for cols in range(len(row1)):
                physical_parameter = row1[cols]
                # print(physical_parameter)
                # planetary_data.setdefault(row[0],{})[physical_parameter] = float(row[cols])
                # print(row[0], physical_parameter, row[cols+1])
                planet_data_dict.setdefault(row[0], {})[physical_parameter] = float(
                    row[cols + 1]
                )
    return planet_data_dict


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


class Venus(Planet):
    """Object containing parameters for Venus."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 0.72  # AU
        self.eccentricity = 0.0067
        self.period = 224.7 / 365.26  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
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


class Mars(Planet):
    """Object containing parameters for Mars."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 1.52  # AU
        self.eccentricity = 0.093
        self.period = 1.88  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()


class Jupiter(Planet):
    """Object containing parameters for Jupiter."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 5.20  # AU
        self.eccentricity = 0.048
        self.period = 11.86  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()


class Saturn(Planet):
    """Object containing parameters for Saturn."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 9.54  # AU
        self.eccentricity = 0.056
        self.period = 29.46  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()


class Uranus(Planet):
    """Object containing parameters for Uranus."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 19.19  # AU
        self.eccentricity = 0.046
        self.period = 84.01  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()


class Neptune(Planet):
    """Object containing parameters for Neptune."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 30.06  # AU
        self.eccentricity = 0.0113
        self.period = 164.79  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()


class Pluto(Planet):
    """Object containing parameters for Pluto."""

    def __init__(self):
        """Physical parameters measured in AU and years."""
        self.mass = None
        self.semimajor = 39.53  # AU
        self.eccentricity = 0.248
        self.period = 248.59  # years
        self.perihelion = self.semimajor * (1 - self.eccentricity)

        self.x0 = self.perihelion + 0j
        self.v0 = self.velocity_at_perihelion()
