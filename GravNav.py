from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math

numpy.seterr(all="print")

###################################################
## Define universal constant of gravitation and
## spacecraft mass (kg), and earth radius (m), along
## with solar constant and speed of light
###################################################
G = numpy.double(6.67384e-20)
mass_spacecraft = numpy.double(100.)
Rearth = numpy.double(6.371e3)
Rsun = numpy.double(6.95800e5)
Rmoon=numpy.double(1.7374e3)
W = numpy.double(1.362e9) #W/sq kilometer
c = numpy.double(2.99792458e5)

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
## to model drag for other bodies, make these fields
## of the planet class and specify them for each
## planet.
###################################################
po = numpy.double(1.01325e2)
To = numpy.double(288.15)
g = numpy.double(9.80665e-3)
L = numpy.double(6.5)
R = numpy.double(8.31447e0)

###################################################
## Define some properties of the spacecraft (these
## are important for modeling drag).
##
## Cd = drag coefficient (assume spherical spacecraft)
## A = effective area (assume 1 sq. meter)
##
## NOTE: Put these in a spacecraft class, makes more
## sense that way.
###################################################
Cd = numpy.double(.47) #sphere
A = numpy.double(1.e-6) #sq m

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
	def __init__(self, mass=0, eci_x=None, eci_y=None, eci_z=None, J2=1.7555e10\
		,x = Symbol('x'), y = Symbol('y'), z = Symbol('z'), potential=None\
		,visualization_potential=None, M = 2.89644e-2, drag_x = None\
		,drag_y=None,drag_z=None, xdot = Symbol('xdot')\
		,ydot = Symbol('ydot'),zdot = Symbol('zdot')):
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

		self.potential = -(((G*self.mass)/sqrt((self.x-self.eci_x)**2. +\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.))\
		+(self.J2*(1./(sqrt(((self.x-self.eci_x)**2.)+((self.y-self.eci_y)**2.)\
			+((self.z-self.eci_z)**2.))**5.))*(1./2.)*\
		(3.*((self.z-self.eci_z)**2.)-(((self.x-self.eci_x)**2.) + \
			(self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.))))

		self.visualization_potential = -(((G*self.mass)/sqrt((self.x-\
			self.eci_x)**2.+(self.y-self.eci_y)**2.))+\
		(self.J2*(1./(sqrt(((self.x-self.eci_x)**2.)+\
			((self.y-self.eci_y)**2.))**5.))*(1./2.)*\
		(-(((self.x-self.eci_x)**2.)+ (self.y-self.eci_y)**2.))))

		#####################################################################
		## NOTE: The drag equations below are true for Earth atmosphere. If 
		## you want to model other planetary atmospheres, make the constants
		## fields for the planet class.
		#####################################################################

		self.drag_x = (.5*(((po*((1.-((L*(sqrt((self.x-self.eci_x)**2.+\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-Rearth))/To))**((g*\
		self.M)/(R*L))))*M)/(R*(To-(L*(sqrt((self.x-self.eci_x)**2.+\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-Rearth)))))*\
		(self.xdot**2.)*Cd*A)

		self.drag_y = (.5*(((po*((1.-((L*(sqrt((self.x-self.eci_x)**2.+\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-Rearth))/To))**((g*\
		self.M)/(R*L))))*M)/(R*(To-(L*(sqrt((self.x-self.eci_x)**2.+\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-Rearth)))))*\
		(self.ydot**2.)*Cd*A)

		self.drag_z = (.5*(((po*((1.-((L*(sqrt((self.x-self.eci_x)**2.+\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-Rearth))/To))**((g*\
		self.M)/(R*L))))*M)/(R*(To-(L*(sqrt((self.x-self.eci_x)**2.+\
		 (self.y-self.eci_y)**2. + (self.z-self.eci_z)**2.)-Rearth)))))*\
		(self.zdot**2.)*Cd*A)

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


#######################################################
## Create solarSystem class, which contains objects
## of the Planet class. For each planet that is added
## to the solar system, the landscape is updated to 
## include that planet's contribution to the potential.
#######################################################

class solarSystem(object):
	def __init__(self, planets = [], landscape = 0, visualization_landscape = 0\
		,vis = None, x = Symbol('x'), y = Symbol('y'), z = Symbol('z')\
		,acc_x = 0, acc_y = 0, acc_z = 0, ext_x = 0, ext_y = 0, ext_z = 0\
		,solarx = 0, solary = 0, solarz = 0,sun_eci_x = 9.583978400242774e7\
		,sun_eci_y = -1.027137783167994e8, sun_eci_z = -4.452859933622356e7):
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
		self.sun_eci_x = sun_eci_x
		self.sun_eci_y = sun_eci_y
		self.sun_eci_z = sun_eci_z

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

		solar_x = self.sun_eci_x - self.x
		solar_y = self.sun_eci_y - self.y
		solar_z = self.sun_eci_z - self.z
		solar_dist = sqrt(solar_x**2. +solar_y**2. +solar_z**2.)

		solar_magnitude = (W/(c*((solar_x)**2. + (solar_y)**2. +\
		(solar_z)**2)*(1./149597870.700)))*A # max angle

		# self.solarx = solar_magnitude*(sqrt(solar_x**2 + solar_y**2)\
		# /solar_dist)*(solar_x/(solar_dist*(sqrt(solar_x**2 + solar_y **2)\
		# /solar_dist)))

		# self.solary = solar_magnitude*(sqrt(solar_x**2 + solar_y**2)\
		# /solar_dist)*(solar_y/(solar_dist*(sqrt(solar_x**2 + solar_y**2)\
		# /solar_dist)))

		self.solarx = solar_magnitude*(solar_x/solar_dist)
		self.solary = solar_magnitude*(solar_y/solar_dist)

		self.solarz = solar_magnitude*(solar_z/solar_dist)

		self.ext_x+=self.solarx
		self.ext_y+=self.solary
		self.ext_z+=self.solarz


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
			self.ext_x/mass_spacecraft)
		self.acc_y += (-diff(self.landscape,self.y) - \
			self.ext_y/mass_spacecraft)
		self.acc_z += (-diff(self.landscape,self.z) - \
			self.ext_z/mass_spacecraft)

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


class Spacecraft(object):
	def __init__(self, xtemp = 0, xdottemp = 0, ytemp = 0, ydottemp = 0\
		, ztemp = 0, zdottemp= 0, xarray = [], yarray = [], zarray = []\
		, xdotarray = [], ydotarray = [], zdotarray = [], xdotn = 0\
		, ydotn = 0, zdotn = 0):
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

	def leapFrog(self,solarsystem,xinit,yinit,zinit,xdotinit,ydotinit,zdotinit,dt,numsteps, x = Symbol('x')\
		, y = Symbol('y'), z = Symbol('z'), xdot = Symbol('xdot')\
		, ydot = Symbol('ydot'), zdot = Symbol('zdot')):
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


		acceleration_x = lambdify((x,y,z,xdot), solarsystem.acc_x\
			, modules = 'numpy')
		acceleration_y = lambdify((x,y,z,ydot), solarsystem.acc_y\
			, modules = 'numpy')
		acceleration_z = lambdify((x,y,z,zdot), solarsystem.acc_z\
			, modules = 'numpy')


		for k in range(numsteps):

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


def lunarCubeSat():
	Earth = Planet(mass=5.972e24,eci_x=0,eci_y=0,eci_z=0,J2=1.7555e10)
	# Moon = Planet(mass=7.34767309e22, eci_x=-1.304659527132289e5,\
	#  eci_y=-3.637250531186822e5, eci_z=-1.221315894487088e5, J2=0., M=0.)

	## THIS IS A FAKE MOON. HEAVY, AND CLOSE
	Moon = Planet(mass=7.34767309e24, eci_x=-2*15115.40312811781-Rearth,\
		eci_y=-2*23668.97680091111-Rearth, eci_z=2*2241.504923500546+Rearth, J2=0., M=0.)

	solarsystem = solarSystem()
	solarsystem.addPlanet(Earth)
	cubesat = Spacecraft()

	fig = plt.figure()
	ax = fig.gca(projection='3d')

	#solarsystem.addPlanet(Moon)

	cubesat.leapFrog(solarsystem, -15015.40312811781-Rearth, -23568.97680091111-Rearth,\
		2241.504923500546+Rearth, -.4855378922082240, -5.048763191594085,\
		-.8799883161637991, 1., 50000)

	ax.plot(numpy.array(cubesat.zarray),\
		numpy.array(cubesat.yarray),\
		numpy.array(cubesat.xarray),\
		label="CubeSat Position, with Moon")

	u = numpy.linspace(0, numpy.pi, 30)
	v = numpy.linspace(0, 2 * numpy.pi, 30)

	xpt = numpy.outer(numpy.sqrt(Rearth)*numpy.sin(u), numpy.sqrt(Rearth)*numpy.sin(v))
	ypt = numpy.outer(numpy.sqrt(Rearth)*numpy.sin(u), numpy.sqrt(Rearth)*numpy.cos(v))
	zpt = numpy.outer(numpy.sqrt(Rearth)*numpy.cos(u), numpy.sqrt(Rearth)*numpy.ones_like(v))


	ax.plot_wireframe(xpt, ypt, zpt,label="Earth")

	xptm = numpy.outer(numpy.sqrt(Rmoon)*numpy.sin(u), numpy.sqrt(Rmoon)*numpy.sin(v))-1.304659527132289e5
	yptm = numpy.outer(numpy.sqrt(Rmoon)*numpy.sin(u), numpy.sqrt(Rmoon)*numpy.cos(v))-3.637250531186822e5
	zptm = numpy.outer(numpy.sqrt(Rmoon)*numpy.cos(u), numpy.sqrt(Rmoon)*numpy.ones_like(v))-1.221315894487088e5

	ax.plot_wireframe(xptm, yptm, zptm,label="Moon")


	plt.title("Passive CubeSat Trajectory")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

def runUnitTest():
	Earth = Planet(mass=5.972e24,eci_x=0,eci_y=0,eci_z=0,J2=1.7555e10)
	Earth_vis = Planet(mass=5.972e50,eci_x=0,eci_y=0,eci_z=0,J2=1.7555e10)
	solarsystem = solarSystem()
	solarsystem.addPlanet(Earth)
	viking = Spacecraft()

	## Planet Potential Landscape
	Earth_vis.showPotential()

	## Falling from height along x axis
	viking.leapFrog(solarsystem,Rearth+.100,0,0,0,0,0,.25,18)
	plt.plot((numpy.array(viking.xarray)-Rearth)*1000,label="position")
	plt.plot(numpy.array(viking.xdotarray)*1000, label="velocity")
	plt.title("x-Position, Falling from 100m Height above Earth's Surface")
	plt.xlabel("time (4.5 sec)")
	plt.ylabel("x")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Falling from height along y axis
	viking.leapFrog(solarsystem,0,Rearth+.100,0,0,0,0,.25,18)
	plt.plot((numpy.array(viking.yarray)-Rearth)*1000,label="position")
	plt.plot(numpy.array(viking.ydotarray)*1000, label="velocity")
	plt.title("y-Position, Falling from 100m Height above Earth's Surface")
	plt.xlabel("time (4.5 sec)")
	plt.ylabel("y")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Falling from height along z axis
	viking.leapFrog(solarsystem,0,0,Rearth+.100,0,0,0,.25,18)
	plt.plot((numpy.array(viking.zarray)-Rearth)*1000,label="position")
	plt.plot(numpy.array(viking.zdotarray)*1000, label="velocity")
	plt.title("z-Position, Falling from 100m Height above Earth's Surface")
	plt.xlabel("time (4.5 sec)")
	plt.ylabel("z")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in xy plane
	viking.leapFrog(solarsystem,Rearth+.50,0,0,0,.500,0,.25,18)
	plt.plot((numpy.array(viking.yarray))*1000,(numpy.array(viking.xarray)-Rearth)*1000\
		,label="ballistic position")
	plt.xlabel("y position")
	plt.ylabel("x position")
	plt.title("Ballistic Trajectory in xy Plane")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in xz plane
	viking.leapFrog(solarsystem,Rearth+.50,0,0,0,0,.500,.25,18)
	plt.plot(numpy.array(viking.zarray)*1000,(numpy.array(viking.xarray)-Rearth)*1000\
		,label="ballistic position")
	plt.xlabel("z position")
	plt.ylabel("x position")
	plt.title("Ballistic Trajectory in xz Plane")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in yz plane
	viking.leapFrog(solarsystem,0,Rearth+.50,0,0,0,.500,.25,18)
	plt.plot(numpy.array(viking.zarray)*1000,(numpy.array(viking.yarray)-Rearth)*1000\
		,label="ballistic position")
	plt.xlabel("z position")
	plt.ylabel("y position")
	plt.title("Ballistic Trajectory in yz Plane")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in xyz space
	viking.leapFrog(solarsystem,Rearth+.50,0,0,0,.500,.500,.25,18)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray)*1000,numpy.array(viking.yarray)*1000\
		,(numpy.array(viking.xarray)-Rearth)*1000,label="ballistic position")
	plt.title("Ballistic Trajectory in 3d Space")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Nodal Regression
	viking.leapFrog(solarsystem, (Rearth), 0., 0., 0., (0),(10),100,2650)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Nodal Regression (short-term)")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Nodal and Apsidal
	viking.leapFrog(solarsystem, (Rearth), 0., 0., 0., (10),(0),200,15000)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Long-Term Orbit around Earth-Like Planet, Nodal and Apsidal Precession")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Just-Less-Than-Escape velocity
	viking.leapFrog(solarsystem, (Rearth), 0., 0., 0., (0),(11.1),100,16500)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Just-Less-Than-Escape Velocity (11.1 km/s)")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Escape velocity
	viking.leapFrog(solarsystem, (Rearth), 0., 0., 0., (0),(11.2),200,16500)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Escape Velocity (11.2 km/s)")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Add a gravitating body
	anotherEarth=Planet(mass=5.972e24,eci_x=0,eci_y=25000,eci_z=20000,J2=1.7555e10)
	solarsystem.addPlanet(anotherEarth)

	## Two Gravitating Bodies
	viking.leapFrog(solarsystem, (Rearth), 0., 0., 0., (10),(0),10,20000)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Trajectory under the influence of Two Earth-Like Planets, 30,000 km Separation")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()


	# ## Solar Pressure
	# sun = solarSystem()
	# viking.leapFrog(solarsystem, solarsystem.sun_eci_x + Rsun\
	# 	, solarsystem.sun_eci_y + Rsun, solarsystem.sun_eci_z + Rsun\
	# 	,0,0,0, 10000, 100000)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot((numpy.array(viking.zdotarray))*1000,\
	# 	(numpy.array(viking.ydotarray))*1000\
	# 	,(numpy.array(viking.xdotarray))*1000\
	# 	,label="solar pressure")
	# plt.title("Spacecraft Velocity from Solar Pressure, Initialized on Solar Surface")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()


	# ## Air damping with altitude - y
	# altitude = []
	# altitude.extend([400])
	# # for i in range(100):
	# # 	altitude.extend([i*10])
	# altitude = numpy.array(altitude)+Rearth
	# for j in altitude:
	# 	viking.leapFrog(solarsystem, j, 0, 0, 0, 5.00, 0, 30, 864000)
	# 	plt.plot(numpy.array(viking.ydotarray)*1000)
	# plt.title("Air Damping with Altitude - y, Surface to 1000 km")
	# plt.xlabel("Time (0-50 seconds)")
	# plt.ylabel("Velocity through Medium (m/s)")
	# plt.show()

	# ## Air damping with altitude - z
	# altitude = []
	# for i in range(100):
	# 	altitude.extend([i*10])
	# altitude = numpy.array(altitude)+Rearth
	# for j in altitude:
	# 	viking.leapFrog(solarsystem, j, 0, 0, 0, 0, 1.00, .25, 200)
	# 	plt.plot(numpy.array(viking.zdotarray)*1000)
	# plt.title("Air Damping with Altitude - z, Surface to 1000 km")
	# plt.xlabel("Time (0-50 seconds)")
	# plt.ylabel("Velocity through Medium (m/s)")
	# plt.show()

	# ## Air damping with altitude - x
	# altitude = []
	# for i in range(100):
	# 	altitude.extend([i*10])
	# altitude = numpy.array(altitude)+Rearth
	# for j in altitude:
	# 	viking.leapFrog(solarsystem, 0, j, 0, 1.00, 0, 0, .25, 200)
	# 	plt.plot(numpy.array(viking.xdotarray)*1000)
	# plt.title("Air Damping with Altitude - x, Surface to 1000 km")
	# plt.xlabel("Time (0-50 seconds)")
	# plt.ylabel("Velocity through Medium (m/s)")
	# plt.show()
	# altitude = []




#lunarCubeSat()
runUnitTest()

## Coming steps:
##	1. Use leapfrog with damping to numerically integrate acc_x,y,z
##	2. Make a sexy plot of the true orbit
##	3. You will be left with an array of positions at each time. Look
##	   up existing technology and use the true position to create a
##	   measured acceleration
##	4. Generate estimated position from measured acceleration
##	5. Go back and Kalman filter measurements
##	6. Add more planets.
##
##	The state will include xarray[], yarray[], zarray[], xdotarray[]
##	ydotarray[], zdotarray[], along with xtemp, ytemp, ztemp, xdottemp
##	ydottemp, and zdottemp


	
