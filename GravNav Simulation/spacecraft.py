from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math
import config
##################################################################
## Create the Spacecraft class ###################################
##################################################################

###################################################
## Define some properties of the spacecraft (these
## are important for modeling drag).
##
## Cd = drag coefficient (assume spherical spacecraft)
## A = effective area (assume 1 sq. meter)
## mass = spacecraft mass (kg)
##
###################################################

class Spacecraft(object):
	def __init__(self, xtemp = 0, xdottemp = 0, ytemp = 0, ydottemp = 0\
		, ztemp = 0, zdottemp= 0, xarray = [], yarray = [], zarray = []\
		, xdotarray = [], ydotarray = [], zdotarray = [], xdotn = 0\
		, ydotn = 0, zdotn = 0, mass = 100., Cd = 0.47, A = numpy.double(1.e-6)):
		self.xtemp = xtemp
		self.xdottemp = xdottemp
		self.ytemp = ytemp
		self.ydottemp = ydottemp
		self.ztemp = ztemp
		self.zdottemp = zdottemp
		self.xarray = xarray
		self.yarray = yarray
		self.zarray = zarray
		self.xdotarray = xdotarray
		self.ydotarray = ydotarray
		self.zdotarray = zdotarray
		self.xdotn = xdotn
		self.ydotn = ydotn
		self.zdotn = zdotn
		self.mass = mass
		self.Cd = Cd
		self.A = A

	def leapFrog(self,solarsystem,xinit,yinit,zinit,xdotinit,ydotinit,zdotinit,dt,numsteps, x = Symbol('x')\
		, y = Symbol('y'), z = Symbol('z'), xdot = Symbol('xdot')\
		, ydot = Symbol('ydot'), zdot = Symbol('zdot'), resetFlag = True\
		, propagationFlag = config.propagationFlag):
		###########################################################
		## Initializations ########################################
		###########################################################
		self.xtemp = 0
		self.ytemp = 0
		self.ztemp = 0
		self.xdottemp = 0
		self.ydottemp = 0
		self.zdottemp = 0
		self.xdotn = 0
		self.ydotn = 0
		self.zdotn = 0

		if resetFlag == True:
			self.xarray = []
			self.yarray = []
			self.zarray = []
			self.xdotarray = []
			self.ydotarray = []
			self.zdotarray = []

			self.xarray.extend([xinit])
			self.xarray.extend([xinit])
			self.yarray.extend([yinit])
			self.yarray.extend([yinit])
			self.zarray.extend([zinit])
			self.zarray.extend([zinit])
			self.xdotarray.extend([xdotinit])
			self.ydotarray.extend([ydotinit])
			self.zdotarray.extend([zdotinit])

		else:
			pass


		acceleration_x = lambdify((x,y,z,xdot), solarsystem.acc_x\
			, modules = 'numpy')
		acceleration_y = lambdify((x,y,z,ydot), solarsystem.acc_y\
			, modules = 'numpy')
		acceleration_z = lambdify((x,y,z,zdot), solarsystem.acc_z\
			, modules = 'numpy')


		for k in range(numsteps):

			if propagationFlag == True:
				a = int((k*dt)/60.)
				b = int(((k-1)*dt)/60.)
				if k>0 and (a != b):
					solarsystem.moveForward(a-b)
					acceleration_x = lambdify((x,y,z,xdot), solarsystem.acc_x\
						, modules = 'numpy')
					acceleration_y = lambdify((x,y,z,ydot), solarsystem.acc_y\
						, modules = 'numpy')
					acceleration_z = lambdify((x,y,z,zdot), solarsystem.acc_z\
						, modules = 'numpy')
				else:
					pass

			self.xdottemp = (self.xdotarray[-1] + dt*\
				(acceleration_x(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.xdotarray[-1]+0j)))
			self.ydottemp = (self.ydotarray[-1] + dt*\
				(acceleration_y(self.xarray[-1]+0j\
				,self.yarray[-1]+0j,self.zarray[-1]+0j,self.ydotarray[-1]+0j)))
			self.zdottemp = (self.zdotarray[-1] + dt*\
				(acceleration_z(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.zdotarray[-1]+0j)))



			self.xtemp = (self.xarray[-1] + dt*self.xdottemp)
			self.ytemp = (self.yarray[-1] + dt*self.ydottemp)
			self.ztemp = (self.zarray[-1] + dt*self.zdottemp)



			self.xdotn = ((self.xtemp - self.xarray[-2])/(2*dt))
			self.ydotn = ((self.ytemp - self.yarray[-2])/(2*dt))
			self.zdotn = ((self.ztemp - self.zarray[-2])/(2*dt))



			self.xdotarray.extend([(self.xdotarray[-1] + \
				dt*acceleration_x(self.xarray[-1]+0j,self.yarray[-1]+0j,\
					self.zarray[-1]+0j,self.xdotn+0j))])
			self.ydotarray.extend([(self.ydotarray[-1] + \
				dt*acceleration_y(self.xarray[-1]+0j,self.yarray[-1]+0j,\
					self.zarray[-1]+0j,self.ydotn+0j))])
			self.zdotarray.extend([(self.zdotarray[-1] + \
				dt*acceleration_z(self.xarray[-1]+0j,self.yarray[-1]+0j,\
					self.zarray[-1]+0j,self.zdotn+0j))])



			self.xarray.extend([(self.xarray[-1] + dt*self.xdotarray[-1])])
			self.yarray.extend([(self.yarray[-1] + dt*self.ydotarray[-1])])
			self.zarray.extend([(self.zarray[-1] + dt*self.zdotarray[-1])])

			percentage = (float(k)/numsteps)*100.
			print str(percentage)+' percent finished'
			print "x: "+str(self.xarray[-1])
			print "y: "+str(self.yarray[-1])
			print "z: "+str(self.zarray[-1])
			print "\n\n"

