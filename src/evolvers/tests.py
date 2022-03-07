from planets import *

semimajor_axis = {
    "Mercury": 0.3871,
    "Venus": 0.7233,
    "Earth": 1.0000,
    "Mars": 1.5237,
    "Jupiter": 5.2028,
    "Saturn": 9.5388,
    "Uranus": 19.182,
    "Neptune": 30.058,
    "Pluto": 39.518,
    "Halley": 18,
}

Planet = Mercury

semimajor_axis[Planet.__name__], getattr(Planet(), "semimajor")
