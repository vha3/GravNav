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
class Planet(object):
	def __init__(self, spacecraft, configuration, index = None, mass=None\
		, eci_x=None, eci_y=None, eci_z=None,x = Symbol('x'), y = Symbol('y')\
		, z = Symbol('z'), potential=None,visualization_potential=None\
		, drag_x = None,drag_y=None,drag_z=None, xdot = Symbol('xdot')\
		,ydot = Symbol('ydot'),zdot = Symbol('zdot'), po =None,To = None\
		, g = None,L=None, R = None,M = None, J2=None, radius = None\
		,dragflag = None, c=None, W=None, G=None, planetx = None\
		,planety = None, planetz=None):

		self.configuration = configuration

		executable = "planetdata = config."+self.configuration
		exec(executable)
		print planetdata

		nextecutable = "universedata = config.universe"
		exec(nextecutable)

		self.index=config.index

		self.mass = planetdata[3]
		self.eci_x = planetdata[0][self.index]
		self.eci_y = planetdata[1][self.index]
		self.eci_z = planetdata[2][self.index]
		self.J2 = planetdata[5]
		self.x = x
		self.y = y
		self.z = z
		self.M = planetdata[11]
		self.xdot = xdot
		self.ydot = ydot
		self.zdot = zdot
		self.po = planetdata[6]
		self.To = planetdata[7]
		self.g = planetdata[8]
		self.L = planetdata[9]
		self.R = planetdata[10]
		self.radius = planetdata[4]
		self.spacecraft = spacecraft
		self.dragflag = planetdata[12]
		self.c = universedata[4]
		self.W = universedata[5]
		self.G = universedata[6]

		self.potential = (-(((self.G*self.mass)/sqrt((self.x-self.eci_x)**2. +\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.))\
		+(self.J2*(1./(sqrt(((self.x-self.eci_x)**2.)+((self.y-self.eci_y)**2.)\
			+((self.z-self.eci_z)**2.))**5.))*(1./2.)*\
		(3.*((self.z-self.eci_z)**2.)-(((self.x-self.eci_x)**2.) + \
			(self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)))))

		if planetdata[0][self.index]!=0:
			self.planetx = (-self.G*self.mass*planetdata[0][self.index])/(((planetdata[0][self.index])**2. +\
				(planetdata[1][self.index])**2. + (planetdata[2][self.index])**2)**(3/2.))

			self.planety = (-self.G*self.mass*planetdata[1][self.index])/(((planetdata[0][self.index])**2. +\
				(planetdata[1][self.index])**2. + (planetdata[2][self.index])**2)**(3/2.))

			self.planetz = (-self.G*self.mass*planetdata[2][self.index])/(((planetdata[0][self.index])**2. +\
				(planetdata[1][self.index])**2. + (planetdata[2][self.index])**2)**(3/2.))

		else:
			self.planetx = 0
			self.planety = 0
			self.planetz = 0

		self.visualization_potential = -(((self.G*self.mass)/sqrt((self.x-\
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

	def reset(self, spacecraft, configuration, index = None, mass=None\
		, eci_x=None, eci_y=None, eci_z=None,x = Symbol('x'), y = Symbol('y')\
		, z = Symbol('z'), potential=None,visualization_potential=None\
		, drag_x = None,drag_y=None,drag_z=None, xdot = Symbol('xdot')\
		,ydot = Symbol('ydot'),zdot = Symbol('zdot'), po =None,To = None\
		, g = None,L=None, R = None,M = None, J2=None, radius = None\
		,dragflag = None, c=None, W=None, G=None, planetx = None\
		,planety = None, planetz=None):

		self.configuration = configuration

		executable = "planetdata = config."+self.configuration
		exec(executable)

		nextecutable = "universedata = config.universe"
		exec(nextecutable)

		self.index=config.index

		self.mass = planetdata[3]
		self.eci_x = planetdata[0][self.index]
		self.eci_y = planetdata[1][self.index]
		self.eci_z = planetdata[2][self.index]
		self.J2 = planetdata[5]
		self.x = x
		self.y = y
		self.z = z
		self.M = planetdata[11]
		self.xdot = xdot
		self.ydot = ydot
		self.zdot = zdot
		self.po = planetdata[6]
		self.To = planetdata[7]
		self.g = planetdata[8]
		self.L = planetdata[9]
		self.R = planetdata[10]
		self.radius = planetdata[4]
		self.spacecraft = spacecraft
		self.dragflag = planetdata[12]
		self.c = universedata[4]
		self.W = universedata[5]
		self.G = universedata[6]

		self.potential = (-(((self.G*self.mass)/sqrt((self.x-self.eci_x)**2. +\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.))\
		+(self.J2*(1./(sqrt(((self.x-self.eci_x)**2.)+((self.y-self.eci_y)**2.)\
			+((self.z-self.eci_z)**2.))**5.))*(1./2.)*\
		(3.*((self.z-self.eci_z)**2.)-(((self.x-self.eci_x)**2.) + \
			(self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)))))

		if planetdata[0][self.index]!=0:
			self.planetx = (-self.G*self.mass*planetdata[0][self.index])/(((planetdata[0][self.index])**2. +\
				(planetdata[1][self.index])**2. + (planetdata[2][self.index])**2)**(3/2.))

			self.planety = (-self.G*self.mass*planetdata[1][self.index])/(((planetdata[0][self.index])**2. +\
				(planetdata[1][self.index])**2. + (planetdata[2][self.index])**2)**(3/2.))

			self.planetz = (-self.G*self.mass*planetdata[2][self.index])/(((planetdata[0][self.index])**2. +\
				(planetdata[1][self.index])**2. + (planetdata[2][self.index])**2)**(3/2.))

		else:
			self.planetx = 0
			self.planety = 0
			self.planetz = 0

		self.visualization_potential = -(((self.G*self.mass)/sqrt((self.x-\
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

