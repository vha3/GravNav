from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math

####################################
## NOTE: Using units of km, NOT m ##
####################################

# Universal constant of graviatation
G = numpy.double(6.67384e-20)

# Solar constant
W=numpy.double(1.362e9)

# Speed of light
c = numpy.double(2.99792458e5)