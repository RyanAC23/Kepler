import mmf_setup

mmf_setup.set_path()
from kepler.planets import *

# For a proper testfile, this should have measurements of all of the values kepler.planets calculates.

# planets = {'planet' : {'mass', 'semimajor', 'eccentricity', 'period', 'perihelion'}}
planets = {
    #    "Mercury" : {'mass': 0.05528 , 'semimajor': 0.3871, 'e': 0.2056, 'period': 0.2408, 'p': },
    #    "Venus"   : { , 'semimajor' : 0.7233, , , },
    #    "Earth"   : { , 'semimajor' : 1.0000, , , },
    #    "Mars"    : { , 'semimajor' : 1.5237, , , },
    #    "Jupiter" : { , 'semimajor' : 5.2028, , , },
    #    "Saturn"  : { , 'semimajor' : 9.5388, , , },
    #    "Uranus"  : { , 'semimajor' : 19.182, , , },
    #    "Neptune" : { , 'semimajor' : 30.058, , , },
    #    "Pluto"   : { , 'semimajor' : 39.518, , , },
    #    "Halley"  : { , 'semimajor' : 18, , , },
}


# -------------------------------------
# OLD
# -------------------------------------
# semimajor_axis = {
#     "Mercury": 0.3871,
#     "Venus": 0.7233,
#     "Earth": 1.0000,
#     "Mars": 1.5237,
#     "Jupiter": 5.2028,
#     "Saturn": 9.5388,
#     "Uranus": 19.182,
#     "Neptune": 30.058,
#     "Pluto": 39.518,
#     "Halley": 18,
# }

# semimajor_axis[Planet.__name__], getattr(Planet(), "semimajor")
