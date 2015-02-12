import config
from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math
###################################################
## Create Planet class. Contains planet mass and
## coordinates in ECI (J2000). From this, the class
## generates an expression for the gravitational
## potential contribution (symbolic) in terms of
## x, y, and z in ECI J2000.
##
## The visualization potential is the planet's
## potential projected onto the xy plane, for
## plotting.
##
## TODO: dynamically update planet position using
##		 ephemerides. check if spacecraft behind
##		 planet.
##
###################################################

###################################################
## Define constants necessary for modeling Earth
## atmospheric drag
##
## po = sea level standard atmospheric pressure
## To = sea level standard temperature
## g = Earth-surface gravitational acceleration
## L = temperature lapse rate
## R = Ideal gas constant
## M = molar mass of dry air (parameter for planet)
##
## NOTE: These values are Earth-specific. If you wish
## to model drag for other bodies, change the values
## from their defaults
###################################################

class Planet(object):
	def __init__(self, spacecraft, mass=0, eci_x=None, eci_y=None, eci_z=None\
		,x = Symbol('x'), y = Symbol('y'), z = Symbol('z'), potential=None\
		,visualization_potential=None, drag_x = None\
		,drag_y=None,drag_z=None, xdot = Symbol('xdot')\
		,ydot = Symbol('ydot'),zdot = Symbol('zdot'), po =numpy.double(1.01325e2)\
		,To = numpy.double(288.15), g = numpy.double(9.80665e-3)\
		,L=numpy.double(6.5), R = numpy.double(8.31447e0)\
		,M = 2.89644e-2, J2=1.7555e10, radius = numpy.double(6.371e3)\
		,dragflag = False):
		self.mass = mass
		self.eci_x = eci_x
		self.eci_y = eci_y
		self.eci_z = eci_z
		self.J2 = J2
		self.x = x
		self.y = y
		self.z = z
		self.M = M
		self.xdot = xdot
		self.ydot = ydot
		self.zdot = zdot
		self.po = po
		self.To = To
		self.g = g
		self.L = L
		self.R = R
		self.radius = radius
		self.spacecraft = spacecraft
		self.dragflag = dragflag

		self.potential = -(((config.G*self.mass)/sqrt((self.x-self.eci_x)**2. +\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.))\
		+(self.J2*(1./(sqrt(((self.x-self.eci_x)**2.)+((self.y-self.eci_y)**2.)\
			+((self.z-self.eci_z)**2.))**5.))*(1./2.)*\
		(3.*((self.z-self.eci_z)**2.)-(((self.x-self.eci_x)**2.) + \
			(self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.))))

		self.visualization_potential = -(((config.G*self.mass)/sqrt((self.x-\
			self.eci_x)**2.+(self.y-self.eci_y)**2.))+\
		(self.J2*(1./(sqrt(((self.x-self.eci_x)**2.)+\
			((self.y-self.eci_y)**2.))**5.))*(1./2.)*\
		(-(((self.x-self.eci_x)**2.)+ (self.y-self.eci_y)**2.))))

		#####################################################################
		## NOTE: The drag equations below are true for Earth atmosphere. If 
		## you want to model other planetary atmospheres, make the constants
		## fields for the planet class.
		#####################################################################
		if self.dragflag == True:
			self.drag_x = (.5*(((self.po*((1.-((self.L*(sqrt((self.x-self.eci_x)**2.+\
			 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-self.radius))/self.To))**((self.g*\
			self.M)/(self.R*self.L))))*M)/(self.R*(self.To-(self.L*(sqrt((self.x-self.eci_x)**2.+\
			 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-self.radius)))))*\
			(self.xdot**2.)*self.spacecraft.Cd*self.spacecraft.A)

			self.drag_y = (.5*(((self.po*((1.-((self.L*(sqrt((self.x-self.eci_x)**2.+\
			 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-self.radius))/self.To))**((self.g*\
			self.M)/(self.R*self.L))))*M)/(self.R*(self.To-(self.L*(sqrt((self.x-self.eci_x)**2.+\
			 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-self.radius)))))*\
			(self.ydot**2.)*self.spacecraft.Cd*self.spacecraft.A)

			self.drag_z = (.5*(((self.po*((1.-((self.L*(sqrt((self.x-self.eci_x)**2.+\
			 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-self.radius))/self.To))**((self.g*\
			self.M)/(self.R*self.L))))*M)/(self.R*(self.To-(self.L*(sqrt((self.x-self.eci_x)**2.+\
			 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-self.radius)))))*\
			(self.zdot**2.)*self.spacecraft.Cd*self.spacecraft.A)

		else:
			self.drag_x = 0.
			self.drag_y = 0.
			self.drag_z = 0.

	def showPotential(self):
		lowerxlim = self.eci_x - 100.000
		upperxlim = self.eci_x + 100.000
		lowerylim = self.eci_y - 100.000
		upperylim = self.eci_y + 100.000
		lowerzlim = -1500
		upperzlim = 0

		self.vis=plot3d_parametric_surface(self.x,self.y,\
			self.visualization_potential,(self.x,-1000,1000)\
			,(self.y,-1000,1000),(self.z,-1000,1000))

