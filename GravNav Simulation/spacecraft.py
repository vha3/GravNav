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
		, xdotarray = [], ydotarray = [], zdotarray = [], xdotdotarray = []\
		, ydotdotarray = [], zdotdotarray = [], xdotn = 0, ydotn = 0\
		, zdotn = 0, mass = None, Cd = None, A = None):
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
		self.xdotdotarray = xdotdotarray
		self.ydotdotarray = ydotdotarray
		self.zdotdotarray = zdotdotarray
		self.xdotn = xdotn
		self.ydotn = ydotn
		self.zdotn = zdotn
		self.mass = config.scmass
		self.Cd = config.Cd
		self.A = config.A

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
			self.xdotdotarray = []
			self.ydotdotarray = []
			self.zdotdotarray = []

			config.index = 0

			self.xarray.extend([xinit])
			#self.xarray.extend([xinit])
			self.yarray.extend([yinit])
			#self.yarray.extend([yinit])
			self.zarray.extend([zinit])
			#self.zarray.extend([zinit])
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


		self.xdotdotarray.extend([acceleration_x(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.xdotarray[-1]+0j)])
		self.ydotdotarray.extend([acceleration_y(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.ydotarray[-1]+0j)])
		self.zdotdotarray.extend([acceleration_z(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.zdotarray[-1]+0j)])
		self.xdotdotarray.extend([acceleration_x(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.xdotarray[-1]+0j)])
		self.ydotdotarray.extend([acceleration_y(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.ydotarray[-1]+0j)])
		self.zdotdotarray.extend([acceleration_z(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.zdotarray[-1]+0j)])



		for k in range(numsteps):

			print "Number of planets: "+str(len(solarsystem.planets))

			if propagationFlag == True:
				a = int((k*dt)/150.)
				b = int(((k-1)*dt)/150.)
				if k>1 and (a != b):
					print "Index number: "+str(config.index)
					solarsystem.moveForward(a-b)
					acceleration_new_x = lambdify((x,y,z,xdot), solarsystem.acc_x\
						, modules = 'numpy')
					acceleration_new_y = lambdify((x,y,z,ydot), solarsystem.acc_y\
						, modules = 'numpy')
					acceleration_new_z = lambdify((x,y,z,zdot), solarsystem.acc_z\
						, modules = 'numpy')


					if numpy.real(acceleration_x(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.xdotarray[-2]+0j)) <\
						numpy.real(acceleration_new_x(self.xarray[-1]+0j\
							,self.yarray[-1]+0j,self.zarray[-1]+0j,self.xdotarray[-1]+0j)):
						new_weightx=1
						old_weightx=0

					else:
						new_weightx=1
						old_weightx=0

					if numpy.real(acceleration_y(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.ydotarray[-2]+0j)) <\
						numpy.real(acceleration_new_y(self.xarray[-1]+0j\
							,self.yarray[-1]+0j,self.zarray[-1]+0j,self.ydotarray[-1]+0j)):
						new_weighty=1
						old_weighty=0

					else:
						new_weighty=1
						old_weighty=0					

					if numpy.real(acceleration_z(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.zdotarray[-2]+0j)) <\
						numpy.real(acceleration_new_z(self.xarray[-1]+0j\
							,self.yarray[-1]+0j,self.zarray[-1]+0j,self.zdotarray[-1]+0j)):
						new_weightz=1
						old_weightz=0

					else:
						new_weightz=1
						old_weightz=0

					# old_weightx=abs(numpy.real(acceleration_x(self.xarray[-2]+0j,self.yarray[-2]+0j\
					# 		,self.zarray[-2]+0j,self.xdotarray[-2]+0j)))/\
					# (abs(numpy.real(acceleration_x(self.xarray[-2]+0j,self.yarray[-2]+0j\
					# 		,self.zarray[-2]+0j,self.xdotarray[-2]+0j))) +\
					# 	abs(numpy.real(acceleration_new_x(self.xarray[-1]+0j\
					# 		,self.yarray[-1]+0j,self.zarray[-1]+0j,self.xdotarray[-1]+0j))))

					# new_weightx = 1.-old_weightx

					# old_weighty=abs(numpy.real(acceleration_y(self.xarray[-2]+0j,self.yarray[-2]+0j\
					# 		,self.zarray[-2]+0j,self.ydotarray[-2]+0j)))/\
					# (abs(numpy.real(acceleration_y(self.xarray[-2]+0j,self.yarray[-2]+0j\
					# 		,self.zarray[-2]+0j,self.ydotarray[-2]+0j))) +\
					# 	abs(numpy.real(acceleration_new_y(self.xarray[-1]+0j\
					# 		,self.yarray[-1]+0j,self.zarray[-1]+0j,self.ydotarray[-1]+0j))))

					# new_weighty = 1.-old_weighty

					# old_weightz=abs(numpy.real(acceleration_z(self.xarray[-2]+0j,self.yarray[-2]+0j\
					# 		,self.zarray[-2]+0j,self.zdotarray[-2]+0j)))/\
					# (abs(numpy.real(acceleration_z(self.xarray[-2]+0j,self.yarray[-2]+0j\
					# 		,self.zarray[-2]+0j,self.zdotarray[-2]+0j))) +\
					# 	abs(numpy.real(acceleration_new_z(self.xarray[-1]+0j\
					# 		,self.yarray[-1]+0j,self.zarray[-1]+0j,self.zdotarray[-1]+0j))))

					# new_weightz = 1.-old_weightz







					self.xdottemp = (self.xdotarray[-2] + dt*\
						((old_weightx*acceleration_x(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.xdotarray[-2]+0j)+\
						new_weightx*acceleration_new_x(self.xarray[-1]+0j\
							,self.yarray[-1]+0j,self.zarray[-1]+0j,self.xdotarray[-1]+0j))))
					self.ydottemp = (self.ydotarray[-2] + dt*\
						((old_weighty*acceleration_y(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.ydotarray[-2]+0j)+\
						new_weighty*acceleration_new_y(self.xarray[-1]+0j\
							,self.yarray[-1]+0j,self.zarray[-1]+0j,self.ydotarray[-1]+0j))))
					self.zdottemp = (self.zdotarray[-2] + dt*\
						((old_weightz*acceleration_z(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.zdotarray[-2]+0j)+\
						new_weightz*acceleration_new_z(self.xarray[-1]+0j\
							,self.yarray[-1]+0j,self.zarray[-1]+0j,self.zdotarray[-1]+0j))))



					self.xtemp = (self.xarray[-2] + dt*self.xdottemp)
					self.ytemp = (self.yarray[-2] + dt*self.ydottemp)
					self.ztemp = (self.zarray[-2] + dt*self.zdottemp)


					if k==1:
						self.xdotn = ((self.xtemp - self.xarray[-2])/(dt))
						self.ydotn = ((self.ytemp - self.yarray[-2])/(dt))
						self.zdotn = ((self.ztemp - self.zarray[-2])/(dt))

					else:
						self.xdotn = ((self.xtemp - self.xarray[-3])/(2*dt))
						self.ydotn = ((self.ytemp - self.yarray[-3])/(2*dt))
						self.zdotn = ((self.ztemp - self.zarray[-3])/(2*dt))



					self.xdotarray[-1]=(self.xdotarray[-2] + \
						dt*((old_weightx*acceleration_x(self.xarray[-2]+0j,self.yarray[-2]+0j,\
							self.zarray[-2]+0j,self.xdotn+0j)+\
						new_weightx*acceleration_new_x(self.xarray[-1]+0j\
						,self.yarray[-1]+0j,self.zarray[-1]+0j,self.xdotn+0j))))

					self.ydotarray[-1]=(self.ydotarray[-2] + \
						dt*((old_weighty*acceleration_y(self.xarray[-2]+0j,self.yarray[-2]+0j,\
							self.zarray[-2]+0j,self.ydotn+0j)+\
						new_weighty*acceleration_new_y(self.xarray[-1]+0j\
						,self.yarray[-1]+0j,self.zarray[-1]+0j,self.ydotn+0j))))

					self.zdotarray[-1]=(self.zdotarray[-2] + \
						dt*((old_weightz*acceleration_z(self.xarray[-2]+0j,self.yarray[-2]+0j,\
							self.zarray[-2]+0j,self.zdotn+0j)+\
						new_weightz*acceleration_new_z(self.xarray[-1]+0j\
						,self.yarray[-1]+0j,self.zarray[-1]+0j,self.zdotn+0j))))



					self.xarray[-1]=(self.xarray[-2] + dt*self.xdotarray[-2])
					self.yarray[-1]=(self.yarray[-2] + dt*self.ydotarray[-2])
					self.zarray[-1]=(self.zarray[-2] + dt*self.zdotarray[-2])


					self.xdotdotarray[-1]=(old_weightx*acceleration_x(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.xdotarray[-2]+0j) + \
					new_weightx*acceleration_new_x(self.xarray[-1]+0j,self.yarray[-1]+0j,self.zarray[-1]+0j\
						,self.xdotarray[-1]+0j))

					self.ydotdotarray[-1]=(old_weighty*acceleration_y(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.ydotarray[-2]+0j) + \
					new_weighty*acceleration_new_y(self.xarray[-1]+0j,self.yarray[-1]+0j,self.zarray[-1]+0j\
						,self.ydotarray[-1]+0j))

					self.zdotdotarray[-1]=(old_weightz*acceleration_z(self.xarray[-2]+0j,self.yarray[-2]+0j\
							,self.zarray[-2]+0j,self.zdotarray[-2]+0j) + \
					new_weightz*acceleration_new_z(self.xarray[-1]+0j,self.yarray[-1]+0j,self.zarray[-1]+0j\
						,self.zdotarray[-1]+0j))


					acceleration_x = lambdify((x,y,z,xdot), solarsystem.acc_x\
						, modules = 'numpy')
					acceleration_y = lambdify((x,y,z,ydot), solarsystem.acc_y\
						, modules = 'numpy')
					acceleration_z = lambdify((x,y,z,zdot), solarsystem.acc_z\
						, modules = 'numpy')

				else:
					pass


######################################################################
########################### DO NOT TOUCH #############################
######################################################################
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


			if k==0:
				self.xdotn = ((self.xtemp - self.xarray[-1])/(dt))
				self.ydotn = ((self.ytemp - self.yarray[-1])/(dt))
				self.zdotn = ((self.ztemp - self.zarray[-1])/(dt))
			else:
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


			self.xdotdotarray.extend([acceleration_x(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.xdotarray[-1]+0j)])
			self.ydotdotarray.extend([acceleration_y(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.ydotarray[-1]+0j)])
			self.zdotdotarray.extend([acceleration_z(self.xarray[-1]+0j,self.yarray[-1]+0j\
					,self.zarray[-1]+0j,self.zdotarray[-1]+0j)])
######################################################################
########################### DO NOT TOUCH #############################
######################################################################

			percentage = (float(k)/numsteps)*100.
			print str(int(percentage))+' percent finished'
			print "x: "+str(self.xarray[-1])
			print "y: "+str(self.yarray[-1])
			print "z: "+str(self.zarray[-1])
			print "\n\n"
