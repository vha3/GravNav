import config
from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math
from copy import deepcopy
#######################################################
## Create solarSystem class, which contains objects
## of the Planet class. For each planet that is added
## to the solar system, the landscape is updated to 
## include that planet's contribution to the potential.
#######################################################
class solarSystem(object):
	def __init__(self, spacecraft, index=None, planets = [], landscape = 0, visualization_landscape = 0\
		,vis = None, x = Symbol('x'), y = Symbol('y'), z = Symbol('z')\
		,acc_x = 0, acc_y = 0, acc_z = 0, ext_x = 0., ext_y = 0., ext_z = 0.\
		,solarx = 0, solary = 0, solarz = 0,sun_eci_x = None\
		,sun_eci_y = None, sun_eci_z = None\
		,solarFlag = None, c = None, W = None, G = None):

		executable = "universedata = config.universe"
		exec(executable)
		self.index = config.index

		self.planets = planets
		self.landscape = landscape
		self.visualization_landscape = visualization_landscape
		self.x = x
		self.y = y
		self.z = z
		self.acc_x = acc_x
		self.acc_y = acc_y
		self.acc_z = acc_z
		self.ext_x = ext_x
		self.ext_y = ext_y
		self.ext_z = ext_z
		self.sun_eci_x = universedata[0][self.index]
		self.sun_eci_y = universedata[1][self.index]
		self.sun_eci_z = universedata[2][self.index]
		self.spacecraft = spacecraft
		self.solarFlag = universedata[3]
		self.c = universedata[4]
		self.W = universedata[5]
		self.G = universedata[6]

		#####################################################################
		## Currently, solar pressure does not check if the spacecraft has a
		## direct line of sight to the sun. It just assumes that the solar
		## radiation is always exerting a force along the vector connecting
		## spacecraft to sun.
		##
		## I also assume that the incident sunlight hits a planar surface, so
		## this is an overestimate. Solar pressure is a feature of the solar
		## system.
		#####################################################################
		if self.solarFlag == True:
			solar_x = self.sun_eci_x - self.x
			solar_y = self.sun_eci_y - self.y
			solar_z = self.sun_eci_z - self.z
			solar_dist = sqrt(solar_x**2. +solar_y**2. +solar_z**2.)

			solar_magnitude = (self.W/(self.c*((solar_x)**2. + (solar_y)**2. +\
			(solar_z)**2)*(1./149597870.700)))*spacecraft.A # max angle

			self.solarx = solar_magnitude*(solar_x/solar_dist)
			self.solary = solar_magnitude*(solar_y/solar_dist)
			self.solarz = solar_magnitude*(solar_z/solar_dist)

		else:
			self.solarx = 0
			self.solary = 0
			self.solarz = 0

		self.ext_x+=self.solarx
		self.ext_y+=self.solary
		self.ext_z+=self.solarz

	def reset(self, spacecraft, index=None, planets = [], landscape = 0, visualization_landscape = 0\
		,vis = None, x = Symbol('x'), y = Symbol('y'), z = Symbol('z')\
		,acc_x = 0, acc_y = 0, acc_z = 0, ext_x = 0., ext_y = 0., ext_z = 0.\
		,solarx = 0, solary = 0, solarz = 0,sun_eci_x = None\
		,sun_eci_y = None, sun_eci_z = None\
		,solarFlag = None, c = None, W = None, G = None):

		executable = "universedata = config.universe"
		exec(executable)
		self.index = config.index

		self.planets = []
		self.landscape = 0
		self.visualization_landscape = 0
		self.x = x
		self.y = y
		self.z = z
		self.acc_x = 0.
		self.acc_y = 0.
		self.acc_z = 0.
		self.ext_x = 0.
		self.ext_y = 0.
		self.ext_z = 0.
		self.sun_eci_x = universedata[0][self.index]
		self.sun_eci_y = universedata[1][self.index]
		self.sun_eci_z = universedata[2][self.index]
		self.spacecraft = spacecraft
		self.solarFlag = universedata[3]
		self.c = universedata[4]
		self.W = universedata[5]
		self.G = universedata[6]

		#####################################################################
		## Currently, solar pressure does not check if the spacecraft has a
		## direct line of sight to the sun. It just assumes that the solar
		## radiation is always exerting a force along the vector connecting
		## spacecraft to sun.
		##
		## I also assume that the incident sunlight hits a planar surface, so
		## this is an overestimate. Solar pressure is a feature of the solar
		## system.
		#####################################################################
		if self.solarFlag == True:
			solar_x = self.sun_eci_x - self.x
			solar_y = self.sun_eci_y - self.y
			solar_z = self.sun_eci_z - self.z
			solar_dist = sqrt(solar_x**2. +solar_y**2. +solar_z**2.)

			solar_magnitude = (self.W/(self.c*((solar_x)**2. + (solar_y)**2. +\
			(solar_z)**2)*(1./149597870.700)))*spacecraft.A # max angle

			self.solarx = solar_magnitude*(solar_x/solar_dist)
			self.solary = solar_magnitude*(solar_y/solar_dist)
			self.solarz = solar_magnitude*(solar_z/solar_dist)

		else:
			self.solarx = 0
			self.solary = 0
			self.solarz = 0

		self.ext_x+=self.solarx
		self.ext_y+=self.solary
		self.ext_z+=self.solarz

	def moveForward(self, amount):
		config.index+=amount
		sc = self.spacecraft
		planetlist = deepcopy(self.planets)
		for i in planetlist:
			conf = i.configuration
			i.reset(sc,conf)
		self.reset(sc)
		for i in planetlist:
			self.addPlanet(i)	

	def addPlanet(self,planet):
		self.planets.extend([planet])
		self.landscape+=planet.potential
		self.visualization_landscape+=planet.visualization_potential

		self.ext_x+=planet.drag_x
		self.ext_y+=planet.drag_y
		self.ext_z+=planet.drag_z

		################################################################
		## NOTE: When you add external forces, they will show them-
		## selves in the below acceleration equations. Numerically
		## integrate these equations to get the motion of the space-
		## craft. Also, if you come up with a way to measure acceleration
		## directly, these form a system of equations that could be used
		## to find the absolute location of the spacecraft.
		################################################################
		self.acc_x += (-diff(self.landscape,self.x) - \
			self.ext_x/self.spacecraft.mass)
		self.acc_y += (-diff(self.landscape,self.y) - \
			self.ext_y/self.spacecraft.mass)
		self.acc_z += (-diff(self.landscape,self.z) - \
			self.ext_z/self.spacecraft.mass)

	def showPotential(self):
		xlist = []
		ylist = []
		for i in self.planets:
			xlist.extend([i.eci_x])
			ylist.extend([i.eci_y])

		minx = min(xlist) - 100.000
		maxx = max(xlist) + 100.000
		miny = min(ylist) - 100.000
		maxy = max(ylist) + 100.000
		minz = -15.00
		maxz = 0.

		self.vis=plot3d_parametric_surface(self.x,self.y,\
			-self.visualization_landscape,(self.x,-1000,1000)\
			,(self.y,-1000,1000),(self.z,-1400000,0))
