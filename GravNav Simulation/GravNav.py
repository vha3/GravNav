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

def runUnitTestStatic():
	##########################################################################
	############################ THE STATIC UNIVERSE #########################
	##########################################################################
	viking = Spacecraft()
	Earth = Planet(spacecraft = viking, configuration="earth")
	heavyEarth = Planet(spacecraft = viking, configuration="heavyearth")
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
	xdirection = numpy.array(viking.yarray)[:]*numpy.array(viking.zdotarray)[:] - numpy.array(viking.zarray)[:]*numpy.array(viking.ydotarray)[:]
	ydirection = numpy.array(viking.zarray)[:]*numpy.array(viking.xdotarray)[:] - numpy.array(viking.xarray)[:]*numpy.array(viking.zdotarray)[:]
	zdirection = numpy.array(viking.xarray)[:]*numpy.array(viking.ydotarray)[:] - numpy.array(viking.yarray)[:]*numpy.array(viking.xdotarray)[:]
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

	## And the energy
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
	anotherxEarth=Planet(spacecraft=viking, configuration="anotherxearth")#mass=5.972e24,eci_x=0,eci_y=25000,eci_z=20000,J2=1.7555e10)
	solarsystem.addPlanet(anotherxEarth)

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


	viking.leapFrog(solarsystem, 50000,0,0,0,0,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - X")
	plt.legend(loc='upper right')
	plt.show()


	viking.leapFrog(solarsystem, 50000,0,0,0,.5,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - X")
	plt.legend(loc='upper right')
	plt.show()

	viking.leapFrog(solarsystem, 50000,0,0,0,0,0.5,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - X")
	plt.legend(loc='upper right')
	plt.show()

	newsolarsystem = solarSystem(spacecraft=viking)
	anotheryEarth=Planet(spacecraft=viking, configuration="anotheryearth")
	newsolarsystem.addPlanet(Earth)
	newsolarsystem.addPlanet(anotheryEarth)

	viking.leapFrog(newsolarsystem, 0,50000,0,0,0,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - Y")
	plt.legend(loc='upper right')
	plt.show()


	viking.leapFrog(newsolarsystem, 0,50000,0,.5,0,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - Y")
	plt.legend(loc='upper right')
	plt.show()

	viking.leapFrog(newsolarsystem, 0,50000,0,0,0,0.5,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - Y")
	plt.legend(loc='upper right')
	plt.show()

	newnewsolarsystem = solarSystem(spacecraft=viking)
	anotherzEarth=Planet(spacecraft=viking, configuration="anotherzearth")
	newnewsolarsystem.addPlanet(Earth)
	newnewsolarsystem.addPlanet(anotherzEarth)

	viking.leapFrog(newnewsolarsystem, 0,0,50000,0,0,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - Z")
	plt.legend(loc='upper right')
	plt.show()


	viking.leapFrog(newnewsolarsystem, 0,0,50000,.5,0,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - Z")
	plt.legend(loc='upper right')
	plt.show()

	viking.leapFrog(newnewsolarsystem, 0,0,50000,0,0.5,0,60,2000)
	plt.plot(viking.xarray,label='x')
	plt.plot(viking.yarray,label='y')
	plt.plot(viking.zarray,label='z')
	plt.xlabel("time")
	plt.ylabel("position")
	plt.title("Adding Planets Verification - Z")
	plt.legend(loc='upper right')
	plt.show()
	##########################################################################
	############################ THE STATIC UNIVERSE #########################
	##########################################################################




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
	viking1 = Spacecraft()
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
	# 	,75,10612)	
	solarsystem.addPlanet(Sun)
	viking1.leapFrog(solarsystem, xpos, ypos, zpos, xvel, yvel, zvel\
		,75,10612)

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.xarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.zarray),\
	# 	label="Sim Trajectory - Sans Sun")
	ax.plot(numpy.array(viking1.xarray),\
		numpy.array(viking1.yarray),\
		numpy.array(viking1.zarray),\
		label="Sim Trajectory - With Sun")
	# ax.plot(numpy.array(config.stkscx5_2x_sans_sun),\
	# 	numpy.array(config.stkscy5_2x_sans_sun),\
	# 	numpy.array(config.stkscz5_2x_sans_sun),\
	# 	label="STK Trajectory - Sans Sun")
	ax.plot(numpy.array(config.stkscx5_2x_with_sun),\
		numpy.array(config.stkscy5_2x_with_sun),\
		numpy.array(config.stkscz5_2x_with_sun),\
		label="STK Trajectory - With Sun")
	ax.plot(numpy.array(config.stkscx5_2x_all),\
		numpy.array(config.stkscy5_2x_all),\
		numpy.array(config.stkscz5_2x_all),\
		label="STK Trajectory - All Effects")
	ax.plot(config.moonx[:config.index],\
		config.moony[:config.index],config.moonz[:config.index],\
		label="Moon Position")

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


	plt.title("Passive Trajectory")
	ax.legend(loc='upper right',prop={'size':10})
	plt.show()

	# plotx=[]
	# ploty=[]
	# plotz=[]
	plotx1=[]
	ploty1=[]
	plotz1=[]
	for i in range(len(viking1.xarray)):
		if i%4==0:
			# plotx.extend([viking.xarray[i]])
			# ploty.extend([viking.yarray[i]])
			# plotz.extend([viking.zarray[i]])
			plotx1.extend([viking1.xarray[i]])
			ploty1.extend([viking1.yarray[i]])
			plotz1.extend([viking1.zarray[i]])

	# viking.xarray=plotx
	# viking.yarray=ploty
	# viking.zarray=plotz
	viking1.xarray=plotx1
	viking1.yarray=ploty1
	viking1.zarray=plotz1


	# plt.plot(viking.xarray, label="sim sans sun")
	plt.plot(viking1.xarray, label="sim with sun")
	#plt.plot(config.stkscx5_2x_sans_sun,label="stk sans sun")
	plt.plot(config.stkscx5_2x_with_sun,label="stk with sun")
	plt.plot(config.stkscx5_2x_all,label="stk, all effects")
	plt.title('Spacecraft x position')
	plt.xlabel('Time (9 days)')
	plt.ylabel('Spacecraft position (km)')
	plt.legend(loc='upper left')
	plt.show()

	# plt.plot(viking.yarray, label="sim sans sun")
	plt.plot(viking1.yarray, label="sim with sun")
	#plt.plot(config.stkscy5_2x_sans_sun,label="stk sans sun")
	plt.plot(config.stkscy5_2x_with_sun,label="stk with sun")
	plt.plot(config.stkscy5_2x_all,label="stk, all effects")
	plt.xlabel('Time (9 days)')
	plt.ylabel('Spacecraft position (km)')
	plt.title('Spacecraft y position')
	plt.legend(loc='lower left')
	plt.show()

	# plt.plot(viking.zarray, label="sim sans sun")
	plt.plot(viking1.zarray, label="sim with sun")
	#plt.plot(config.stkscz5_2x_sans_sun,label="stk sans sun")
	plt.plot(config.stkscz5_2x_with_sun,label="stk with sun")
	plt.plot(config.stkscz5_2x_all,label="stk, all effects")
	plt.xlabel('Time (9 days)')
	plt.ylabel('Spacecraft position (km)')
	plt.title('Spacecraft z position');
	plt.legend(loc='upper left')
	plt.show()


	# plt.plot(numpy.array(viking.xarray)[:2652]-\
	# 	numpy.array(config.stkscx5_2x_sans_sun)[:2652],\
	# 	 label="error in km - sans sun sim & sans sun STK")
	# plt.plot(numpy.array(viking.xarray)[:2652]-\
	# 	 numpy.array(config.stkscx5_2x_with_sun)[:2652],\
	# 	 label="error in km - sans sun sim & with sun STK")
	# plt.plot(numpy.array(viking1.xarray)[:2652]-\
	# 	numpy.array(config.stkscx5_2x_sans_sun)[:2652],\
	# 	 label="error in km - with sun sim & sans sun STK")
	plt.plot(numpy.array(viking1.xarray)[:2652]-\
		 numpy.array(config.stkscx5_2x_with_sun)[:2652],\
		 label="error in km - with sun sim & with sun STK")
	# plt.plot(numpy.array(viking.xarray)[:2652]-\
	# 	numpy.array(config.stkscx5_2x_all)[:2652],\
	# 	 label="error in km - sans sun sim & all effects STK")
	# plt.plot(numpy.array(viking.xarray)[:2652]-\
	# 	 numpy.array(config.stkscx5_2x_all)[:2652],\
	# 	 label="error in km - sans sun sim & all effects STK")
	# plt.plot(numpy.array(viking1.xarray)[:2652]-\
	# 	numpy.array(config.stkscx5_2x_all)[:2652],\
	# 	 label="error in km - with sun sim & all effects STK")
	plt.plot(numpy.array(viking1.xarray)[:2652]-\
		 numpy.array(config.stkscx5_2x_all)[:2652],\
		 label="error in km - with sun sim & all effects STK")

	plt.title('Accumulated Error between Python and STK 9 days - x')
	plt.xlabel('Time (9 days)')
	plt.ylabel('Error in km')
	plt.legend(loc='lower left')
	plt.show()



	# plt.plot(numpy.array(viking.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_sans_sun)[:2652],\
	# 	 label="error in km - sans sun sim & sans sun STK")
	# plt.plot(numpy.array(viking.yarray)[:2652]-\
	# 	 numpy.array(config.stkscy5_2x_with_sun)[:2652],\
	# 	 label="error in km - sans sun sim & with sun")
	# plt.plot(numpy.array(viking1.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_sans_sun)[:2652],\
	# 	 label="error in km - with sun sim & sans sun STK")
	plt.plot(numpy.array(viking1.yarray)[:2652]-\
		 numpy.array(config.stkscy5_2x_with_sun)[:2652],\
		 label="error in km - with sun sim & with sun")
	# plt.plot(numpy.array(viking.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_all)[:2652],\
	# 	 label="error in km - sans sun sim & all effects STK")
	# plt.plot(numpy.array(viking.yarray)[:2652]-\
	# 	 numpy.array(config.stkscy5_2x_all)[:2652],\
	# 	 label="error in km - sans sun sim & all effects STK")
	# plt.plot(numpy.array(viking1.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_all)[:2652],\
	# 	 label="error in km - with sun sim & all effects STK")
	plt.plot(numpy.array(viking1.yarray)[:2652]-\
		 numpy.array(config.stkscy5_2x_all)[:2652],\
		 label="error in km - with sun sim & all effects STK")
	plt.xlabel('Time (9 days)')
	plt.ylabel('Error in km')
	plt.title('Accumulated Error between Python and STK 9 days - y')
	plt.legend(loc='lower left')
	plt.show()



	# plt.plot(numpy.array(viking.zarray)[:2652]-\
	# 	numpy.array(config.stkscz5_2x_sans_sun)[:2652],\
	# 	 label="error in km - sans sun sim & sans sun STK")
	# plt.plot(numpy.array(viking.zarray)[:2652]-\
	# 	 numpy.array(config.stkscz5_2x_with_sun)[:2652],\
	# 	 label="error in km - sans sun sim & with sun STK")
	# plt.plot(numpy.array(viking1.zarray)[:2652]-\
	# 	numpy.array(config.stkscz5_2x_sans_sun)[:2652],\
	# 	 label="error in km - with sun sim & sans sun STK")
	plt.plot(numpy.array(viking1.zarray)[:2652]-\
		 numpy.array(config.stkscz5_2x_with_sun)[:2652],\
		 label="error in km - with sun sim & with sun STK")
	# plt.plot(numpy.array(viking.zarray)[:2652]-\
	# 	numpy.array(config.stkscz5_2x_all)[:2652],\
	# 	 label="error in km - sans sun sim & all effects STK")
	# plt.plot(numpy.array(viking.zarray)[:2652]-\
	# 	 numpy.array(config.stkscz5_2x_all)[:2652],\
	# 	 label="error in km - sans sun sim & all effects STK")
	# plt.plot(numpy.array(viking1.zarray)[:2652]-\
	# 	numpy.array(config.stkscz5_2x_all)[:2652],\
	# 	 label="error in km - with sun sim & all effects STK")
	plt.plot(numpy.array(viking1.zarray)[:2652]-\
		 numpy.array(config.stkscz5_2x_all)[:2652],\
		 label="error in km - with sun sim & all effects STK")
	plt.xlabel('Time (9 days)')
	plt.ylabel('Error in km')
	plt.title('Accumulated Error between Python and STK 9 days - z')
	plt.legend(loc='lower left')
	plt.show()



	# plt.plot(numpy.sqrt((numpy.array(viking.xarray)[:2652]\
	# 	-numpy.array(config.stkscx5_2x_sans_sun)[:2652])**2 + \
	# (numpy.array(viking.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_sans_sun)[:2652])**2 +\
	# 	 (numpy.array(viking.zarray)[:2652]-\
	# 	 	numpy.array(config.stkscz5_2x_sans_sun)[:2652])**2)\
	# ,label="Total Accumulated Error Sans Sun Sim & Sans Sun STK- km")
	# plt.plot(numpy.sqrt((numpy.array(viking.xarray)[:2652]\
	# 	-numpy.array(config.stkscx5_2x_with_sun)[:2652])**2 + \
	# (numpy.array(viking.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_with_sun)[:2652])**2 +\
	# 	 (numpy.array(viking.zarray)[:2652]-\
	# 	 	numpy.array(config.stkscz5_2x_with_sun)[:2652])**2)\
	# ,label="Total Accumulated Error Sans Sun Sim & With Sun STK- km")
	# plt.plot(numpy.sqrt((numpy.array(viking1.xarray)[:2652]\
	# 	-numpy.array(config.stkscx5_2x_sans_sun)[:2652])**2 + \
	# (numpy.array(viking1.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_sans_sun)[:2652])**2 +\
	# 	 (numpy.array(viking1.zarray)[:2652]-\
	# 	 	numpy.array(config.stkscz5_2x_sans_sun)[:2652])**2)\
	# ,label="Total Accumulated Error With Sun Sim & Sans Sun STK- km")
	plt.plot(numpy.sqrt((numpy.array(viking1.xarray)[:2652]\
		-numpy.array(config.stkscx5_2x_with_sun)[:2652])**2 + \
	(numpy.array(viking1.yarray)[:2652]-\
		numpy.array(config.stkscy5_2x_with_sun)[:2652])**2 +\
		 (numpy.array(viking1.zarray)[:2652]-\
		 	numpy.array(config.stkscz5_2x_with_sun)[:2652])**2)\
	,label="Total Accumulated Error With Sun Sim & With Sun STK- km")
	# plt.plot(numpy.sqrt((numpy.array(viking.xarray)[:2652]\
	# 	-numpy.array(config.stkscx5_2x_all)[:2652])**2 + \
	# (numpy.array(viking.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_all)[:2652])**2 +\
	# 	 (numpy.array(viking.zarray)[:2652]-\
	# 	 	numpy.array(config.stkscz5_2x_all)[:2652])**2)\
	# ,label="Total Accumulated Error Sans Sun Sim & All Effects STK- km")
	# plt.plot(numpy.sqrt((numpy.array(viking.xarray)[:2652]\
	# 	-numpy.array(config.stkscx5_2x_all)[:2652])**2 + \
	# (numpy.array(viking.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_all)[:2652])**2 +\
	# 	 (numpy.array(viking.zarray)[:2652]-\
	# 	 	numpy.array(config.stkscz5_2x_all)[:2652])**2)\
	# ,label="Total Accumulated Error Sans Sun Sim & All Effects STK- km")
	# plt.plot(numpy.sqrt((numpy.array(viking1.xarray)[:2652]\
	# 	-numpy.array(config.stkscx5_2x_all)[:2652])**2 + \
	# (numpy.array(viking1.yarray)[:2652]-\
	# 	numpy.array(config.stkscy5_2x_all)[:2652])**2 +\
	# 	 (numpy.array(viking1.zarray)[:2652]-\
	# 	 	numpy.array(config.stkscz5_2x_all)[:2652])**2)\
	# ,label="Total Accumulated Error With Sun Sim & All Effects STK- km")
	plt.plot(numpy.sqrt((numpy.array(viking1.xarray)[:2652]\
		-numpy.array(config.stkscx5_2x_all)[:2652])**2 + \
	(numpy.array(viking1.yarray)[:2652]-\
		numpy.array(config.stkscy5_2x_all)[:2652])**2 +\
		 (numpy.array(viking1.zarray)[:2652]-\
		 	numpy.array(config.stkscz5_2x_all)[:2652])**2)\
	,label="Total Accumulated Error With Sun Sim & All Effects STK- km")
	plt.title('Total Accumulated Error - 9 days')
	plt.legend(loc='upper left')
	plt.show()



cubeSat()
#runUnitTestStatic()
