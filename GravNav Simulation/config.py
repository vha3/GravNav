from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math

#Position of moon every 5 min for 20 hr
stkmoonx5, stkmoony5, stkmoonz5 = numpy.loadtxt("STK_5min_moon.txt", unpack=True)
#position of spacecraft every 5 min for 20 hr
stkscx5, stkscy5, stkscz5 = numpy.loadtxt("STK_5min_sc.txt", unpack=True)
#position of moon every 5 min for 40 hr
stkmoonx5_2x, stkmoony5_2x, stkmoonz5_2x = numpy.loadtxt("STK_moon_5min_2x.txt", unpack=True)
#position of sun every 5 min for 40 hr
stksunx5_2x, stksuny5_2x, stksunz5_2x = numpy.loadtxt("STK_sun_5min_2x.txt", unpack=True)
#position of spacecraft, sans sun drag and pressure, for 40 hr
stkscx5_2x_sans_sun, stkscy5_2x_sans_sun, stkscz5_2x_sans_sun = numpy.loadtxt("STK_sat_no_drag_press_sun.txt", unpack=True)
#position of spacecraft, sans drag and pressure, for 40 hr
stkscx5_2x_with_sun, stkscy5_2x_with_sun, stkscz5_2x_with_sun = numpy.loadtxt("STK_sat_no_drag_press.txt", unpack=True)
#position of moon every 1 min for 20 hr
stkmoonx1, stkmoony1, stkmoonz1 = numpy.loadtxt("STK_1min_moon.txt", unpack=True)
#position of spacecraft every 1 min for 20 hr (all effects)
stkscx1, stkscy1, stkscz1 = numpy.loadtxt("STK_1min_sc.txt", unpack=True)
#position of spacecraft every 5 min for 40 hour (all effects)
stkscx5_2x_all, stkscy5_2x_all, stkscz5_2x_all = numpy.loadtxt("STK_satellite_2x_with_drag_sun_pressure.txt", unpack=True)

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
# t=Symbol('t')
# moonxeqn = (1.9e-25)*(t**5)-(2.2e-19)*(t**4)-(8.7e-13)*(t**3)\
# +(7.2e-7)*(t**2)+.78*t-2.3e5
# moonyeqn = (4.5e-26)*(t**5)-(5.6e-19)*(t**4)+(5.1e-13)*(t**3)\
# +(9.5e-7)*(t**2)-.53*t-3.1e5
# moonzeqn = (3.4e-27)*(t**5)-(1.8e-19)*(t**4)+(2.4e-13)*(t**3)\
# +(2.9e-7)*(t**2)-.24*t-9.8e4
# xeq = lambdify((t),moonxeqn)
# yeq = lambdify((t),moonyeqn)
# zeq = lambdify((t),moonzeqn)

# moonx = []
# moony = []
# moonz = []
# for i in range(2654):
# 	moonx.extend([xeq(i*300)])
# 	moony.extend([yeq(i*300)])
# 	moonz.extend([zeq(i*300)])

moonx, moony, moonz = numpy.loadtxt("STK_moon_5min_2x.txt", unpack=True) #km



moonxinterp=[]
moonyinterp=[]
moonzinterp=[]
for i in range(len(moonx)-1):
  moonxinterp.extend([(moonx[i]+moonx[i+1])/2])
  moonyinterp.extend([(moony[i]+moony[i+1])/2])
  moonzinterp.extend([(moonz[i]+moonz[i+1])/2])

newmoonx=[]
newmoony=[]
newmoonz=[]

for i in range(len(moonx)-1):
  newmoonx.extend([moonx[i]])
  newmoonx.extend([moonxinterp[i]])
  newmoony.extend([moony[i]])
  newmoony.extend([moonyinterp[i]])
  newmoonz.extend([moonz[i]])
  newmoonz.extend([moonzinterp[i]])



moonmass = 7.34767309e22 #kg
moonradius = 1737.53 #km
moonJ2 = 0.
moonpo = 0.
moonTo = 0.
moong = 0.
moonL = 0.
moonR = 0.
moonM = 0.
moondragFlag = False

moon = [newmoonx, newmoony, newmoonz, moonmass, moonradius,\
 moonJ2, moonpo, moonTo, moong, moonL, moonR, moonM,\
  moondragFlag]

########################################################
## The Sun #############################################
########################################################
# sunxeqn = (-2.2e-13)*(t**3)+(3.5e-7)*(t**2)+30*t-1.6e7
# sunyeqn = (-1.6e-20)*(t**4)+(2.9e-14)*(t**3)+(2.8e-6)*(t**2)\
# -3*t-(1.3e8)
# sunzeqn = (-6.6e-21)*(t**4)+(1.3e-14)*(t**3)+(1.2e-6)*(t**2)\
# -1.3*t-(5.8e7)
# sunxeq = lambdify((t),sunxeqn)
# sunyeq = lambdify((t),sunyeqn)
# sunzeq = lambdify((t),sunzeqn)

# sunx = []
# suny = []
# sunz = []
# for i in range(2654):
#   sunx.extend([sunxeq(i*300)])
#   suny.extend([sunyeq(i*300)])
#   sunz.extend([sunzeq(i*300)])

sunx, suny, sunz = numpy.loadtxt("STK_sun_5min_2x.txt", unpack=True)


sunxinterp=[]
sunyinterp=[]
sunzinterp=[]
for i in range(len(sunx)-1):
  sunxinterp.extend([(sunx[i]+sunx[i+1])/2])
  sunyinterp.extend([(suny[i]+suny[i+1])/2])
  sunzinterp.extend([(sunz[i]+sunz[i+1])/2])

newsunx=[]
newsuny=[]
newsunz=[]

for i in range(len(sunx)-1):
  newsunx.extend([sunx[i]])
  newsunx.extend([sunxinterp[i]])
  newsuny.extend([suny[i]])
  newsuny.extend([sunyinterp[i]])
  newsunz.extend([sunz[i]])
  newsunz.extend([sunzinterp[i]])



sunmass = 1.988544e30
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

sun = [newsunx, newsuny, newsunz, sunmass, sunradius,\
 sunJ2, sunpo, sunTo, sung, sunL, sunR, sunM,\
  sundragFlag]

########################################################
## The Earth #############################################
########################################################
earthx, earthy, earthz = numpy.array(numpy.zeros(len(moon[0]))),\
numpy.array(numpy.zeros(len(moon[0]))),numpy.array(numpy.zeros(len(moon[0])))
earthmass = 5.97219e24
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
## Earth that propagates through space
########################################################
earthpx = []
earthpy = []
earthpz = []
for i in range(len(earthx)):
	earthpx.extend([10*i])
	earthpy.extend([10*i])
	earthpz.extend([10*i])

earth2 = [earthpx, earthpy, earthpz, earthmass, earthradius,\
 earthJ2, earthpo, earthTo, earthg, earthL, earthR, earthM,\
  earthdragFlag]



########################################################
## Heavy Earth #########################################
########################################################
heavyearthx, heavyearthy, heavyearthz = numpy.array(numpy.zeros(len(moon[0]))),\
numpy.array(numpy.zeros(len(moon[0]))),numpy.array(numpy.zeros(len(moon[0])))
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
anotherxearthx, anotherxearthy, anotherxearthz = numpy.array(numpy.zeros(len(moon[0]))),\
numpy.array(numpy.zeros(len(moon[0]))),numpy.array(numpy.zeros(len(moon[0])))
anotherxearthx = anotherxearthx + 100000
anotherxearthmass = 5.97219e24
anotherxearthradius = numpy.double(6.371e3)
anotherxearthJ2 = numpy.double(1.7555e10)
anotherxearthpo = numpy.double(1.01325e2)
anotherxearthTo = numpy.double(288.15)
anotherxearthg = numpy.double(9.80665e-3)
anotherxearthL = numpy.double(6.5)
anotherxearthR = numpy.double(8.31447e0)
anotherxearthM = 2.89644e-2
anotherxearthdragFlag = False

anotherxearth = [anotherxearthx, anotherxearthy, anotherxearthz, anotherxearthmass, anotherxearthradius,\
 anotherxearthJ2, anotherxearthpo, anotherxearthTo, anotherxearthg, anotherxearthL, anotherxearthR, anotherxearthM,\
  anotherxearthdragFlag]

anotheryearthx, anotheryearthy, anotheryearthz = numpy.array(numpy.zeros(len(moon[0]))),\
numpy.array(numpy.zeros(len(moon[0]))),numpy.array(numpy.zeros(len(moon[0])))
anotheryearthy = anotheryearthy + 100000
anotheryearthmass = 5.97219e24
anotheryearthradius = numpy.double(6.371e3)
anotheryearthJ2 = numpy.double(1.7555e10)
anotheryearthpo = numpy.double(1.01325e2)
anotheryearthTo = numpy.double(288.15)
anotheryearthg = numpy.double(9.80665e-3)
anotheryearthL = numpy.double(6.5)
anotheryearthR = numpy.double(8.31447e0)
anotheryearthM = 2.89644e-2
anotheryearthdragFlag = False

anotheryearth = [anotheryearthx, anotheryearthy, anotheryearthz, anotheryearthmass, anotheryearthradius,\
 anotheryearthJ2, anotheryearthpo, anotheryearthTo, anotheryearthg, anotheryearthL, anotheryearthR, anotheryearthM,\
  anotheryearthdragFlag]


anotherzearthx, anotherzearthy, anotherzearthz = numpy.array(numpy.zeros(len(moon[0]))),\
numpy.array(numpy.zeros(len(moon[0]))),numpy.array(numpy.zeros(len(moon[0])))
anotherzearthz = anotherzearthz + 100000
anotherzearthmass = 5.97219e24
anotherzearthradius = numpy.double(6.371e3)
anotherzearthJ2 = numpy.double(1.7555e10)
anotherzearthpo = numpy.double(1.01325e2)
anotherzearthTo = numpy.double(288.15)
anotherzearthg = numpy.double(9.80665e-3)
anotherzearthL = numpy.double(6.5)
anotherzearthR = numpy.double(8.31447e0)
anotherzearthM = 2.89644e-2
anotherzearthdragFlag = False

anotherzearth = [anotherzearthx, anotherzearthy, anotherzearthz, anotherzearthmass, anotherzearthradius,\
 anotherzearthJ2, anotherzearthpo, anotherzearthTo, anotherzearthg, anotherzearthL, anotherzearthR, anotherzearthM,\
  anotherzearthdragFlag]


#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
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
# km to AU conversion
kmAU = 149597870.700

universe = [sun[0],sun[1],sun[2],solarFlag,c,W,G,kmAU]

####################################
## Spacecraft Parameters ###########
####################################
scmass = 100
Cd = 2.2
A = numpy.double(1.e-6)

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