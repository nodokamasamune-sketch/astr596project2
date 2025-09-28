This generates a stellar cluster of population N using the Kroupa IMF and Plummer sphere models. Users are then able to select a specific integration method, duration of simulation, and timestep, to simulate the dynamics of this cluster. It is also possible to set your own initial positions, velocities, and masses.      

Import classes and functions from the .py scripts in the src folder. Additionally import numpy as np and matplotlib.pyplot as plt.       

To generate masses, input the number of bodies, N,  to generate masses for in the 'imf' function. This will generate an array of shape N of masses. To generate initial positions and velocities for these masses, feed the array of masses into the 'plummer' function. This will generate two arrays of shape N x 3: initial positions and initial velocities, both in cartesian coordinates.     

To simulate the dynamics, call the ODE class. The inputs are:
Method: str, choose: 'euler', 'rk2', 'rk4', 'leapfrog' to select integration approximation method
N: int, number of bodies in simulation
initial_positions: np.array of initial positions in cartesian coordinates
initial_velocities: np.array of initial velocities in cartesian coordinates
Dimensions: int, choose 2D or 3D
duration: float, time in years
dt: float, timestep in years   

This will output the positions, velocities, kinetic energy, gravitational potential energy, total energy, and virial ratio for each object at each integration step.

Various plotting and testing functions are defined in the test folder.
