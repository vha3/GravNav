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
## Planetary positions, minute-by-minute from deployment
## until 1 week after deployment
#########################################################
## po = sea level standard atmospheric pressure
## To = sea level standard temperature
## g = Earth-surface gravitational acceleration
## L = temperature lapse rate
## R = Ideal gas constant
## M = molar mass of dry air (parameter for planet)
#########################################################
##
#########################################################
## The Moon #############################################
#########################################################
moonx, moony, moonz = numpy.loadtxt("moon.txt", unpack=True) #km
moonmass = 734.9e20 #kg
moonradius = 1737.53 #km
moonJ2 = 0.
moonpo = 0.
moonTo = 0.
moong = 0.
moonL = 0.
moonR = 0.
moonM = 0.
moondragFlag = False

moon = [moonx, moony, moonz, moonmass, moonradius,\
 moonJ2, moonpo, moonTo, moong, moonL, moonR, moonM,\
  moondragFlag]

########################################################
## The Sun #############################################
########################################################
sunx, suny, sunz = numpy.loadtxt("sun.txt", unpack=True)
sunmass = 1.988e30
sunradius = 6.963e5
sunJ2 = 0.
sunpo = 0.
sunTo = 0.
sung = 0.
sunL = 0.
sunR = 0.
sunM = 0.
sundragFlag = False
sun = [sunx, suny, sunz]

sun = [sunx, suny, sunz, sunmass, sunradius,\
 sunJ2, sunpo, sunTo, sung, sunL, sunR, sunM,\
  sundragFlag]

########################################################
## The Earth #############################################
########################################################
earthx, earthy, earthz = numpy.array(numpy.zeros(len(moonx))),\
numpy.array(numpy.zeros(len(moonx))),numpy.array(numpy.zeros(len(moonx)))
earthmass = 5.972e24
earthradius = numpy.double(6.371e3)
earthJ2 = numpy.double(1.7555e10)
earthpo = numpy.double(1.01325e2)
earthTo = numpy.double(288.15)
earthg = numpy.double(9.80665e-3)
earthL = numpy.double(6.5)
earthR = numpy.double(8.31447e0)
earthM = 2.89644e-2
earthdragFlag = False

earth = [earthx, earthy, earthz, earthmass, earthradius,\
 earthJ2, earthpo, earthTo, earthg, earthL, earthR, earthM,\
  earthdragFlag]


########################################################
## Heavy Earth #########################################
########################################################
heavyearthx, heavyearthy, heavyearthz = numpy.array(numpy.zeros(len(moonx))),\
numpy.array(numpy.zeros(len(moonx))),numpy.array(numpy.zeros(len(moonx)))
heavyearthmass = 5.972e50
heavyearthradius = numpy.double(6.371e3)
heavyearthJ2 = numpy.double(1.7555e10)
heavyearthpo = numpy.double(1.01325e2)
heavyearthTo = numpy.double(288.15)
heavyearthg = numpy.double(9.80665e-3)
heavyearthL = numpy.double(6.5)
heavyearthR = numpy.double(8.31447e0)
heavyearthM = 2.89644e-2
heavyearthdragFlag = False

heavyearth = [heavyearthx, heavyearthy, heavyearthz, heavyearthmass, heavyearthradius,\
 heavyearthJ2, heavyearthpo, heavyearthTo, heavyearthg, heavyearthL, heavyearthR, heavyearthM,\
  heavyearthdragFlag]  


########################################################
## Another Earth #######################################
########################################################
anotherearthx, anotherearthy, anotherearthz = numpy.array(numpy.zeros(len(moonx))),\
numpy.array(numpy.zeros(len(moonx))),numpy.array(numpy.zeros(len(moonx)))
anotherearthy = anotherearthy + 25000
anotherearthz = anotherearthz + 20000
anotherearthmass = 5.972e24
anotherearthradius = numpy.double(6.371e3)
anotherearthJ2 = numpy.double(1.7555e10)
anotherearthpo = numpy.double(1.01325e2)
anotherearthTo = numpy.double(288.15)
anotherearthg = numpy.double(9.80665e-3)
anotherearthL = numpy.double(6.5)
anotherearthR = numpy.double(8.31447e0)
anotherearthM = 2.89644e-2
anotherearthdragFlag = False

anotherearth = [anotherearthx, anotherearthy, anotherearthz, anotherearthmass, anotherearthradius,\
 anotherearthJ2, anotherearthpo, anotherearthTo, anotherearthg, anotherearthL, anotherearthR, anotherearthM,\
  anotherearthdragFlag]



####################################
## Universal constants #############
## and solar system constants ######
####################################
# Universal constant of graviatation
G = numpy.double(6.67384e-20)
# Speed of light
c = numpy.double(2.99792458e5)
# solar constant
W=numpy.double(1.362e9)
# solar flag turns on/off solar pressure
solarFlag = False

universe = [sunx,suny,sunz,solarFlag,c,W,G]


####################################
## Simulation Parameters ###########
####################################
## Propagate Planets? ##############
####################################
propagationFlag = True
## Particular moment in time. Index 0
## corresponds to Dec 15, 2017 14:56
## and each successive index adds 60
## seconds
index = 0






