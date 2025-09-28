#uniform masses
#random positions in sphere of radius 100 AU
#zero initial velocities
#test with N =10
def random_n_plotting():
    '''
    Plotting  10 randomly generated body interaction. Produces :
        3 panel figure showing energy components vs. time, relative energy error vs. time, and virial ratio q over 100 years for each method : euler,  RK2, RK4, and leapfrog. 

    Parameters
    -------
    None: N = 10; other parameters already preset.
    
    Returns
    -------
    3 panel figure showing energy components vs. time, relative energy error vs. time, and virial ratio q over time for each method : euler,  RK2, RK4, and leapfrog. 
    '''
    import sys
    import os
    module_dir = os.path.abspath('../src') #directory to source

    sys.path.append(module_dir)
    
    from ODE import ODE
    import numpy as np
    import matplotlib.pyplot as plt

    #defining random location generator
    def star_generator_random(N):
        #spherical coordinates
        thetas = np.random.randint(0, 90, N)
        phis = np.random.randint(0, 180, N) 
        radii = np.random.randint(0, 100, N) 

        #cartesian coordinates
        x = radii*np.sin(phis)*np.cos(thetas)
        y = radii*np.sin(phis)*np.sin(thetas)
        z = radii*np.cos(phis)

        positions = np.array([x, y, z])

        initial_positions = np.empty([N, 3])
        for i in range(N):
            coord = positions[0:, i]
            initial_positions[i,:] = coord
        
        initial_velocities = np.zeros([N, 3]) #initial velocities = 0
        masses = np.ones([N, 1]) #masses = 1 solar mass
        return initial_positions, initial_velocities, masses

    ini_pos, ini_vel, masses = star_generator_random(10) #N = 10

    #integrating N=10 gravitational interactions using RK4 and leapfrog
    rk4 = ODE('rk4', 10, ini_pos, ini_vel, 3, 100, 0.01, masses, R = 100) #radius = 10
    leapfrog = ODE('leapfrog', 10, ini_pos, ini_vel, 3, 100, 0.01, masses, R = 100)

    #time
    time = np.linspace(0, 100, 10000)
    
    #relative energy error = (E_tot - E_tot_i) / np.abs(E_tot_i)
    rk4_error = (rk4.E_tot - rk4.E_tot[1]) / np.abs(rk4.E_tot[1])
    leapfrog_error = (leapfrog.E_tot - leapfrog.E_tot[1]) / np.abs(leapfrog.E_tot[1])

    #3 panels per method

    #define figure, axes, and figsize
    fig, axes = plt.subplots(3, 1, figsize=(10,14))

    #first plot: kinetic energy, total energy,and potential energy vs. time
    axes[0].plot(time, rk4.E_tot, linestyle='-.', color='purple', label='RK4 total energy')
    axes[0].plot(time, leapfrog.E_tot, linestyle='-.', color='red', label='Leapfrog total energy')
    axes[0].plot(time, rk4.kE, linestyle='--', color='blue', label='RK4 kinetic energy')
    axes[0].plot(time, leapfrog.kE, linestyle='--', color='green', label='Leapfrog kinetic energy')
    axes[0].plot(time, rk4.U_g, linestyle=':', color='orange', label='RK4 gravitational potential energy')
    axes[0].plot(time, leapfrog.U_g, linestyle=':', color='black', label='Leapfrog gravitational potential energy')
    axes[0].set_title('Total Energy Change over Time')
    axes[0].set_ylabel('Energy (M_solar AU^2/yr^2)')
    axes[0].set_xlabel('Time (yrs)')
    axes[0].grid(True)
    axes[0].legend()    

    axes[1].plot(time, rk4_error, color='purple', linestyle='-.',label='RK4')
    axes[1].plot(time, leapfrog_error, color='red', linestyle=':', label='Leapfrog')
    #axes[1].set_yscale('log')
    axes[1].set_title('Relative Energy Error Change over Time')
    axes[1].set_ylabel('Relative Energy Error')
    axes[1].set_xlabel('Time (yrs)')
    axes[1].grid(True)
    axes[1].legend()    

    axes[2].plot(time, rk4.virial, color='purple', linestyle='-.', label='RK4')
    axes[2].plot(time, leapfrog.virial, color='red', linestyle=':', label='Leapfrog')
    axes[2].set_title('Virial Change over Time')
    axes[2].set_ylabel('Virial Ratio')
    axes[2].set_yscale('log')
    axes[2].set_xlabel('Time (yrs)')
    axes[2].grid(True)
    axes[2].legend()    

    fig.suptitle('Comparing RK4 and Leapfrog Approximation Methods for 10 uniformly randomly generated bodies')
    plt.tight_layout()
    plt.savefig('../outputs/figures/Comparing RK4 and Leapfrog Approximation Methods for 10 uniformly randomly generated bodies')
    plt.show()