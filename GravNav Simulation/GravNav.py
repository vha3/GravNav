from sympy import sqrt, Symbol, diff
from sympy.utilities.lambdify import lambdify
from sympy.plotting import plot3d_parametric_surface
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import math
from spacecraft import Spacecraft
from solarsystem import solarSystem
from planet import Planet
import config

numpy.seterr(all="print")

def runUnitTest():
	viking = Spacecraft()
	Earth = Planet(spacecraft = viking, configuration="earth")# mass=5.972e24,eci_x=0,eci_y=0,eci_z=0)
	heavyEarth = Planet(spacecraft = viking, configuration="heavyearth")# mass=5.972e50,eci_x=0,eci_y=0,eci_z=0)
	solarsystem = solarSystem(spacecraft=viking)
	solarsystem.addPlanet(Earth)


	## Planet Potential Landscape
	heavyEarth.showPotential()

	## Falling from height along x axis
	viking.leapFrog(solarsystem,Earth.radius+.100,0,0,0,0,0,.25,18)
	plt.plot((numpy.array(viking.xarray)-Earth.radius)*1000,label="position")
	plt.plot(numpy.array(viking.xdotarray)*1000, label="velocity")
	plt.title("x-Position, Falling from 100m Height above Earth's Surface")
	plt.xlabel("time (4.5 sec)")
	plt.ylabel("x")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Falling from height along y axis
	viking.leapFrog(solarsystem,0,Earth.radius+.100,0,0,0,0,.25,18)
	plt.plot((numpy.array(viking.yarray)-Earth.radius)*1000,label="position")
	plt.plot(numpy.array(viking.ydotarray)*1000, label="velocity")
	plt.title("y-Position, Falling from 100m Height above Earth's Surface")
	plt.xlabel("time (4.5 sec)")
	plt.ylabel("y")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Falling from height along z axis
	viking.leapFrog(solarsystem,0,0,Earth.radius+.100,0,0,0,.25,18)
	plt.plot((numpy.array(viking.zarray)-Earth.radius)*1000,label="position")
	plt.plot(numpy.array(viking.zdotarray)*1000, label="velocity")
	plt.title("z-Position, Falling from 100m Height above Earth's Surface")
	plt.xlabel("time (4.5 sec)")
	plt.ylabel("z")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in xy plane
	viking.leapFrog(solarsystem,Earth.radius+.50,0,0,0,.500,0,.25,18)
	plt.plot((numpy.array(viking.yarray))*1000,(numpy.array(viking.xarray)-Earth.radius)*1000\
		,label="ballistic position")
	plt.xlabel("y position")
	plt.ylabel("x position")
	plt.title("Ballistic Trajectory in xy Plane")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in xz plane
	viking.leapFrog(solarsystem,Earth.radius+.50,0,0,0,0,.500,.25,18)
	plt.plot(numpy.array(viking.zarray)*1000,(numpy.array(viking.xarray)-Earth.radius)*1000\
		,label="ballistic position")
	plt.xlabel("z position")
	plt.ylabel("x position")
	plt.title("Ballistic Trajectory in xz Plane")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in yz plane
	viking.leapFrog(solarsystem,0,Earth.radius+.50,0,0,0,.500,.25,18)
	plt.plot(numpy.array(viking.zarray)*1000,(numpy.array(viking.yarray)-Earth.radius)*1000\
		,label="ballistic position")
	plt.xlabel("z position")
	plt.ylabel("y position")
	plt.title("Ballistic Trajectory in yz Plane")
	plt.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Ballistic trajectory in xyz space
	viking.leapFrog(solarsystem,Earth.radius+.50,0,0,0,.500,.500,.25,18)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray)*1000,numpy.array(viking.yarray)*1000\
		,(numpy.array(viking.xarray)-Earth.radius)*1000,label="ballistic position")
	plt.title("Ballistic Trajectory in 3d Space")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Nodal Regression
	viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (5.5921),(5.5921),60,270)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Nodal Regression (short-term)")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	## Plot direction of anguluar momentum vector
	viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (5.5921),(5.5921),60,73571)
	xdirection = numpy.array(viking.yarray)[:-1]*numpy.array(viking.zdotarray)[:] - numpy.array(viking.zarray)[:-1]*numpy.array(viking.ydotarray)[:]
	ydirection = numpy.array(viking.zarray)[:-1]*numpy.array(viking.xdotarray)[:] - numpy.array(viking.xarray)[:-1]*numpy.array(viking.zdotarray)[:]
	zdirection = numpy.array(viking.xarray)[:-1]*numpy.array(viking.ydotarray)[:] - numpy.array(viking.yarray)[:-1]*numpy.array(viking.xdotarray)[:]
	for i in range(len(xdirection)):
		xdirection[i] = xdirection[i]/(numpy.sqrt(xdirection[i]**2 + ydirection[i]**2 + zdirection[i]**2))
		ydirection[i] = ydirection[i]/(numpy.sqrt(xdirection[i]**2 + ydirection[i]**2 + zdirection[i]**2))
		zdirection[i] = zdirection[i]/(numpy.sqrt(xdirection[i]**2 + ydirection[i]**2 + zdirection[i]**2))
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(xdirection,ydirection,zdirection,label="Direction of angular momentum vector")
	plt.title("Precession of Angular Momentum Vector for Circular Orbit at Earth Radius, Theoretical Time Required for Full 1 Period")
	ax.legend(loc='upper right')
	plt.show()


	## Nodal and Apsidal
	viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (7),(7),200,15000)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Long-Term Orbit around Earth-Like Planet, Nodal and Apsidal Precession")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()


	## Get plot of energy
	viking.leapFrog(solarsystem, Earth.radius, 0., 0., 0., 7, 7, 1, 40000)
	x = Symbol('x')
	y = Symbol('y')
	z = Symbol('z')
	potential = lambdify((x,y,z), Earth.potential\
			, modules = 'numpy')
	potential_energy = viking.mass*(potential(numpy.array(viking.xarray),numpy.array(viking.yarray),numpy.array(viking.zarray)))
	kinetic_energy = .5*viking.mass*(numpy.array(viking.xdotarray)**2 + numpy.array(viking.ydotarray)**2 + numpy.array(viking.zdotarray)**2)
	total = potential_energy[:40000]+kinetic_energy[:40000]
	plt.plot(potential_energy,label="Potential Energy")
	plt.plot(kinetic_energy,label="Kinetic Energy")
	plt.plot(total,label="Total Energy")
	plt.title("Conservation of Energy (Verification of Numeric Integrator)")
	plt.xlabel("Time (1.5 million timesteps)")
	plt.ylabel("Energy")
	plt.legend()
	plt.show()
	print len(potential_energy)
	print len(kinetic_energy)
	print "Initial total energy: "+str(total[0])+"\n"
	print "Final total energy: "+str(total[-1])+"\n"

	## Just-Less-Than-Escape velocity
	viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (0),(11.1),100,16500)
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
	viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (0),(11.2),200,16500)
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
	anotherEarth=Planet(spacecraft=viking, configuration="anotherearth")#mass=5.972e24,eci_x=0,eci_y=25000,eci_z=20000,J2=1.7555e10)
	solarsystem.addPlanet(anotherEarth)

	## Two Gravitating Bodies
	viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (12),(0),10,20000)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(numpy.array(viking.zarray),\
		numpy.array(viking.yarray),\
		numpy.array(viking.xarray),\
		label="Spacecraft Position")
	plt.title("Trajectory under the influence of Two Earth-Like Planets, 30,000 km Separation")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()


	## Conservation of energy, two bodies
	x = Symbol('x')
	y = Symbol('y')
	z = Symbol('z')
	potential = lambdify((x,y,z), solarsystem.landscape\
			, modules = 'numpy')
	potential_energy = viking.mass*(potential(numpy.array(viking.xarray),numpy.array(viking.yarray),numpy.array(viking.zarray)))
	kinetic_energy = .5*viking.mass*(numpy.array(viking.xdotarray)**2 + numpy.array(viking.ydotarray)**2 + numpy.array(viking.zdotarray)**2)
	total = potential_energy[:20000]+kinetic_energy[:20000]
	plt.plot(potential_energy,label="Potential Energy")
	plt.plot(kinetic_energy,label="Kinetic Energy")
	plt.plot(total,label="Total Energy")
	plt.title("Conservation of Energy, 2 Gravitating Bodies (Verification of Numeric Integrator)")
	plt.xlabel("Time (200,000 timesteps)")
	plt.ylabel("Energy")
	plt.legend()
	plt.show()
	print len(potential_energy)
	print len(kinetic_energy)
	print "Initial total energy: "+str(total[0])+"\n"
	print "Final total energy: "+str(total[-1])+"\n"



	# ## Solar Pressure
	# sun = solarSystem()
	# viking.leapFrog(solarsystem, solarsystem.sun_eci_x + 6.95800e5\
	# 	, solarsystem.sun_eci_y + 6.95800e5, solarsystem.sun_eci_z + 6.95800e5\
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
	# altitude = numpy.array(altitude)+Earth.radius
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
	# altitude = numpy.array(altitude)+Earth.radius
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
	# altitude = numpy.array(altitude)+Earth.radius
	# for j in altitude:
	# 	viking.leapFrog(solarsystem, 0, j, 0, 1.00, 0, 0, .25, 200)
	# 	plt.plot(numpy.array(viking.xdotarray)*1000)
	# plt.title("Air Damping with Altitude - x, Surface to 1000 km")
	# plt.xlabel("Time (0-50 seconds)")
	# plt.ylabel("Velocity through Medium (m/s)")
	# plt.show()
	# altitude = []

def cubeSat():
	viking = Spacecraft()
	Earth = Planet(spacecraft = viking, configuration="earth")
	Moon = Planet(spacecraft = viking, configuration="moon")
	Sun = Planet(spacecraft = viking, configuration="sun")
	solarsystem = solarSystem(spacecraft=viking)
	solarsystem.addPlanet(Earth)
	solarsystem.addPlanet(Moon)
	#solarsystem.addPlanet(Sun)

	## CubeSat Trajectory
	xpos = -1.501540312811781e4
	ypos = -2.356897680091111e4
	zpos =  2.241504923500546e3
	xvel = -4.855378922082240e-1
	yvel = -5.048763191594085e0
	zvel = -8.799883161637991e-1
	# viking.leapFrog(solarsystem, xpos, ypos, zpos, xvel, yvel, zvel\
	# 	,300,242)

	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# ax.plot(config.moonx[:config.index],\
	# 	config.moony[:config.index],config.moonz[:config.index],\
	# 	label="Moon Position")

	# u = numpy.linspace(0, 2 * numpy.pi, 100)
	# v = numpy.linspace(0, numpy.pi, 100)
	# x = Earth.radius * numpy.outer(numpy.cos(u), numpy.sin(v))
	# y = Earth.radius * numpy.outer(numpy.sin(u), numpy.sin(v))
	# z = Earth.radius * numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

	# xm = Moon.radius * numpy.outer(numpy.cos(u), numpy.sin(v))\
	#  + config.moonx[config.index]
	# ym = Moon.radius * numpy.outer(numpy.cos(u), numpy.sin(v))\
	#  + config.moony[config.index]
	# zm = Moon.radius * numpy.outer(numpy.cos(u), numpy.sin(v))\
	#  + config.moonz[config.index]

	# ax.plot_surface(xm,ym,zm, rstride=4, cstride=4,color='r',\
	# 	label="Moon")
	# ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r',\
	#  label="Earth")


	# plt.title("Passive Trajectory")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# plt.plot(viking.xarray, label="sim");plt.plot(config.stkscx5, label="stk"); plt.legend(loc='upper right'); plt.show()
	# plt.plot(viking.yarray, label="sim");plt.plot(config.stkscy5, label="stk"); plt.legend(loc='upper right'); plt.show()
	# plt.plot(viking.zarray, label="sim");plt.plot(config.stkscz5, label="stk"); plt.legend(loc='upper right'); plt.show()

	# plt.plot(config.moonx[:config.index], label="sim");plt.plot(config.stkmoonx5, label="stk"); plt.legend(loc='upper right'); plt.show()
	# plt.plot(config.moony[:config.index], label="sim");plt.plot(config.stkmoony5, label="stk"); plt.legend(loc='upper right'); plt.show()
	# plt.plot(config.moonz[:config.index], label="sim");plt.plot(config.stkmoonz5, label="stk"); plt.legend(loc='upper right'); plt.show()


	viking.leapFrog(solarsystem, xpos, ypos, zpos, xvel, yvel, zvel\
		,30,2402)

	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# ax.plot(config.moonx[:config.index],\
	# 	config.moony[:config.index],config.moonz[:config.index],\
	# 	label="Moon Position")

	# u = numpy.linspace(0, 2 * numpy.pi, 100)
	# v = numpy.linspace(0, numpy.pi, 100)
	# x = Earth.radius * numpy.outer(numpy.cos(u), numpy.sin(v))
	# y = Earth.radius * numpy.outer(numpy.sin(u), numpy.sin(v))
	# z = Earth.radius * numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

	# xm = Moon.radius * numpy.outer(numpy.cos(u), numpy.sin(v))\
	#  + config.moonx[config.index]
	# ym = Moon.radius * numpy.outer(numpy.cos(u), numpy.sin(v))\
	#  + config.moony[config.index]
	# zm = Moon.radius * numpy.outer(numpy.cos(u), numpy.sin(v))\
	#  + config.moonz[config.index]

	# ax.plot_surface(xm,ym,zm, rstride=4, cstride=4,color='r',\
	# 	label="Moon")
	# ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r',\
	#  label="Earth")


	# plt.title("Passive Trajectory")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()

	plotx=[]
	ploty=[]
	plotz=[]
	for i in range(len(viking.xarray)):
		if i%2==0:
			plotx.extend([viking.xarray[i]])
			ploty.extend([viking.yarray[i]])
			plotz.extend([viking.zarray[i]])

	viking.xarray=plotx
	viking.yarray=ploty
	viking.zarray=plotz


	plt.plot(viking.xarray, label="sim");plt.plot(config.stkscx1, label="stk");plt.title('Spacecraft x position'); plt.legend(loc='upper right'); plt.show()
	plt.plot(viking.yarray, label="sim");plt.plot(config.stkscy1, label="stk");plt.title('Spacecraft y position'); plt.legend(loc='upper right'); plt.show()
	plt.plot(viking.zarray, label="sim");plt.plot(config.stkscz1, label="stk");plt.title('Spacecraft z position'); plt.legend(loc='upper right'); plt.show()

	plt.plot(numpy.array(viking.xarray)[:1200]-numpy.array(config.stkscx1)[:1200], label="error in km");plt.title('Accumulated Error between Python and STK 20 hrs - x');plt.legend(loc='upper left'); plt.show()
	plt.plot(numpy.array(viking.yarray)[:1200]-numpy.array(config.stkscy1)[:1200], label="error in km");plt.title('Accumulated Error between Python and STK 20 hrs - y');plt.legend(loc='upper left'); plt.show()
	plt.plot(numpy.array(viking.zarray)[:1200]-numpy.array(config.stkscz1)[:1200], label="error in km");plt.title('Accumulated Error between Python and STK 20 hrs - z');plt.legend(loc='upper left'); plt.show()

	plt.plot(numpy.sqrt((numpy.array(viking.xarray)[:1200]-numpy.array(config.stkscx1)[:1200])**2 + (numpy.array(viking.yarray)[:1200]-numpy.array(config.stkscy1)[:1200])**2 + (numpy.array(viking.zarray)[:1200]-numpy.array(config.stkscz1)[:1200])**2),label="Total Accumulated Error - km");plt.title('Total Accumulated Error - 20 hrs');plt.legend(loc='upper left');plt.show()

	# plt.plot(config.moonx[:config.index], label="sim");plt.plot(config.stkmoonx1, label="stk");plt.title('Moon x Position'); plt.legend(loc='upper right'); plt.show()
	# plt.plot(config.moony[:config.index], label="sim");plt.plot(config.stkmoony1, label="stk");plt.title('Moon y Position'); plt.legend(loc='upper right'); plt.show()
	# plt.plot(config.moonz[:config.index], label="sim");plt.plot(config.stkmoonz1, label="stk");plt.title('Moon z Position'); plt.legend(loc='upper right'); plt.show()

def drifters():

	config.propagationFlag = False


	viking1 = Spacecraft()
	viking2 = Spacecraft()
	Earth = Planet(spacecraft = viking1, configuration="earth")
	solarsystem = solarSystem(spacecraft=viking1)
	solarsystem.addPlanet(Earth)

	viking1.leapFrog(solarsystem, 2*Earth.radius+.001, 0, 0, 0, 0, 0, 1, 300)
	viking2.leapFrog(solarsystem, 2*Earth.radius, 0, 0, 0, 0, 0, 1, 300)

	xdiff = 1000*(numpy.array(viking1.xarray) - numpy.array(viking2.xarray))
	ydiff = 1000*(numpy.array(viking1.yarray) - numpy.array(viking2.yarray))
	zdiff = 1000*(numpy.array(viking1.zarray) - numpy.array(viking2.zarray))

	plt.plot(xdiff,label="x difference")
	plt.plot(ydiff,label="y difference")
	plt.plot(zdiff,label="z difference")

	plt.title("Drift Positions, Dropped from X")
	plt.xlabel("Time (5 min)")
	plt.ylabel("Difference in Position (m)")
	plt.legend(loc='upper right')
	plt.show()

	xdiff = 1000*(numpy.array(viking1.xdotarray) - numpy.array(viking2.xdotarray))
	ydiff = 1000*(numpy.array(viking1.ydotarray) - numpy.array(viking2.ydotarray))
	zdiff = 1000*(numpy.array(viking1.zdotarray) - numpy.array(viking2.zdotarray))

	plt.plot(xdiff,label="x difference")
	plt.plot(ydiff,label="y difference")
	plt.plot(zdiff,label="z difference")

	plt.title("Drift Velocities, Dropped from X")
	plt.xlabel("Time (5 min)")
	plt.ylabel("Difference in Position (m/sec)")
	plt.legend(loc='upper right')
	plt.show()

	xdiff = 1000*(numpy.array(viking1.xdotdotarray) - numpy.array(viking2.xdotdotarray))
	ydiff = 1000*(numpy.array(viking1.ydotdotarray) - numpy.array(viking2.ydotdotarray))
	zdiff = 1000*(numpy.array(viking1.zdotdotarray) - numpy.array(viking2.zdotdotarray))

	plt.plot(xdiff,label="x difference")
	plt.plot(ydiff,label="y difference")
	plt.plot(zdiff,label="z difference")

	plt.title("Relative Acceleration, Dropped from X")
	plt.xlabel("Time (5 min)")
	plt.ylabel("Difference in Acc (m/sec/sec)")
	plt.legend(loc='upper right')
	plt.show()


	viking1.leapFrog(solarsystem, -Earth.radius, 0, 0, 0, -8, 0, 60, 300)
	viking2.leapFrog(solarsystem, Earth.radius, 0, 0, 0, 8, 0, 60, 300)

	xdiff = 1000*(numpy.array(viking1.xarray) - numpy.array(viking2.xarray))
	ydiff = 1000*(numpy.array(viking1.yarray) - numpy.array(viking2.yarray))
	zdiff = 1000*(numpy.array(viking1.zarray) - numpy.array(viking2.zarray))

	plt.plot(xdiff,label="x difference")
	plt.plot(ydiff,label="y difference")
	plt.plot(zdiff,label="z difference")

	plt.title("Drift Positions, Nearly Circular Orbit, Opposite Sides of Earth")
	plt.xlabel("Time (5 min)")
	plt.ylabel("Difference in Position (m)")
	plt.legend(loc='upper right')
	plt.show()

	plt.plot(numpy.sqrt((xdiff**2) + (ydiff**2) + (zdiff**2)));plt.show()


	xdiff = 1000*(numpy.array(viking1.xdotarray) - numpy.array(viking2.xdotarray))
	ydiff = 1000*(numpy.array(viking1.ydotarray) - numpy.array(viking2.ydotarray))
	zdiff = 1000*(numpy.array(viking1.zdotarray) - numpy.array(viking2.zdotarray))

	plt.plot(xdiff,label="x difference")
	plt.plot(ydiff,label="y difference")
	plt.plot(zdiff,label="z difference")

	plt.title("Drift Velocities, Nearly Circular Orbit, Opposite Sides of Earth")
	plt.xlabel("Time (5 min)")
	plt.ylabel("Difference in Position (m/sec)")
	plt.legend(loc='upper right')
	plt.show()

	xdiff = 1000*(numpy.array(viking1.xdotdotarray) - numpy.array(viking2.xdotdotarray))
	ydiff = 1000*(numpy.array(viking1.ydotdotarray) - numpy.array(viking2.ydotdotarray))
	zdiff = 1000*(numpy.array(viking1.zdotdotarray) - numpy.array(viking2.zdotdotarray))

	plt.plot(xdiff,label="x difference")
	plt.plot(ydiff,label="y difference")
	plt.plot(zdiff,label="z difference")

	plt.title("Relative Acceleration, Nearly Circular Orbit, Opposite Sides of Earth")
	plt.xlabel("Time (5 min)")
	plt.ylabel("Difference in Acc (m/sec/sec)")
	plt.legend(loc='upper right')
	plt.show()


cubeSat()
#runUnitTest()
#drifters()

	
