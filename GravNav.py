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

numpy.seterr(all="print")

def runUnitTest():
	viking = Spacecraft()
	Earth = Planet(spacecraft = viking, mass=5.972e24,eci_x=0,eci_y=0,eci_z=0)
	Earth_vis = Planet(spacecraft = viking, mass=5.972e50,eci_x=0,eci_y=0,eci_z=0)
	solarsystem = solarSystem(spacecraft=viking)
	solarsystem.addPlanet(Earth)


	# ## Planet Potential Landscape
	# Earth_vis.showPotential()

	# ## Falling from height along x axis
	# viking.leapFrog(solarsystem,Earth.radius+.100,0,0,0,0,0,.25,18)
	# plt.plot((numpy.array(viking.xarray)-Earth.radius)*1000,label="position")
	# plt.plot(numpy.array(viking.xdotarray)*1000, label="velocity")
	# plt.title("x-Position, Falling from 100m Height above Earth's Surface")
	# plt.xlabel("time (4.5 sec)")
	# plt.ylabel("x")
	# plt.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Falling from height along y axis
	# viking.leapFrog(solarsystem,0,Earth.radius+.100,0,0,0,0,.25,18)
	# plt.plot((numpy.array(viking.yarray)-Earth.radius)*1000,label="position")
	# plt.plot(numpy.array(viking.ydotarray)*1000, label="velocity")
	# plt.title("y-Position, Falling from 100m Height above Earth's Surface")
	# plt.xlabel("time (4.5 sec)")
	# plt.ylabel("y")
	# plt.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Falling from height along z axis
	# viking.leapFrog(solarsystem,0,0,Earth.radius+.100,0,0,0,.25,18)
	# plt.plot((numpy.array(viking.zarray)-Earth.radius)*1000,label="position")
	# plt.plot(numpy.array(viking.zdotarray)*1000, label="velocity")
	# plt.title("z-Position, Falling from 100m Height above Earth's Surface")
	# plt.xlabel("time (4.5 sec)")
	# plt.ylabel("z")
	# plt.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Ballistic trajectory in xy plane
	# viking.leapFrog(solarsystem,Earth.radius+.50,0,0,0,.500,0,.25,18)
	# plt.plot((numpy.array(viking.yarray))*1000,(numpy.array(viking.xarray)-Earth.radius)*1000\
	# 	,label="ballistic position")
	# plt.xlabel("y position")
	# plt.ylabel("x position")
	# plt.title("Ballistic Trajectory in xy Plane")
	# plt.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Ballistic trajectory in xz plane
	# viking.leapFrog(solarsystem,Earth.radius+.50,0,0,0,0,.500,.25,18)
	# plt.plot(numpy.array(viking.zarray)*1000,(numpy.array(viking.xarray)-Earth.radius)*1000\
	# 	,label="ballistic position")
	# plt.xlabel("z position")
	# plt.ylabel("x position")
	# plt.title("Ballistic Trajectory in xz Plane")
	# plt.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Ballistic trajectory in yz plane
	# viking.leapFrog(solarsystem,0,Earth.radius+.50,0,0,0,.500,.25,18)
	# plt.plot(numpy.array(viking.zarray)*1000,(numpy.array(viking.yarray)-Earth.radius)*1000\
	# 	,label="ballistic position")
	# plt.xlabel("z position")
	# plt.ylabel("y position")
	# plt.title("Ballistic Trajectory in yz Plane")
	# plt.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Ballistic trajectory in xyz space
	# viking.leapFrog(solarsystem,Earth.radius+.50,0,0,0,.500,.500,.25,18)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray)*1000,numpy.array(viking.yarray)*1000\
	# 	,(numpy.array(viking.xarray)-Earth.radius)*1000,label="ballistic position")
	# plt.title("Ballistic Trajectory in 3d Space")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Nodal Regression
	# viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (7),(7),100,2650)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# plt.title("Nodal Regression (short-term)")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Nodal and Apsidal
	# viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (7),(7),200,15000)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# plt.title("Long-Term Orbit around Earth-Like Planet, Nodal and Apsidal Precession")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()


	## Get plot of energy
	viking.leapFrog(solarsystem, Earth.radius, 0., 0., 0., 7, 7, 1, 1500000)
	x = Symbol('x')
	y = Symbol('y')
	z = Symbol('z')
	potential = lambdify((x,y,z), Earth.potential\
			, modules = 'numpy')
	potential_energy = viking.mass*(potential(numpy.array(viking.xarray),numpy.array(viking.yarray),numpy.array(viking.zarray)))
	kinetic_energy = .5*viking.mass*(numpy.array(viking.xdotarray)**2 + numpy.array(viking.ydotarray)**2 + numpy.array(viking.zdotarray)**2)
	total = potential_energy[:1500000]+kinetic_energy[:1500000]
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

	# ## Just-Less-Than-Escape velocity
	# viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (0),(11.1),100,16500)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# plt.title("Just-Less-Than-Escape Velocity (11.1 km/s)")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Escape velocity
	# viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (0),(11.2),200,16500)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# plt.title("Escape Velocity (11.2 km/s)")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()

	# ## Add a gravitating body
	# anotherEarth=Planet(spacecraft=viking,mass=5.972e24,eci_x=0,eci_y=25000,eci_z=20000,J2=1.7555e10)
	# solarsystem.addPlanet(anotherEarth)

	# ## Two Gravitating Bodies
	# viking.leapFrog(solarsystem, (Earth.radius), 0., 0., 0., (10),(0),10,20000)
	# fig = plt.figure()
	# ax = fig.gca(projection='3d')
	# ax.plot(numpy.array(viking.zarray),\
	# 	numpy.array(viking.yarray),\
	# 	numpy.array(viking.xarray),\
	# 	label="Spacecraft Position")
	# plt.title("Trajectory under the influence of Two Earth-Like Planets, 30,000 km Separation")
	# ax.legend(loc='upper right',prop={'size':10})
	# plt.show()


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



runUnitTest()


	
