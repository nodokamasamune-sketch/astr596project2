def jupiter_plotting(N, duration, time_step):
    '''
    Plotting 2 body earth, jupiter, sun interaction. Produces 2 x 2 panel comparing trajectories, total energy vs time, relative energy error vs time, and the virial ratio q over time for each method.

    Parameters
    ----------    
    N : float
        Number of bodies to calculate interactions for

    duration : int
        Number of years to simulate.
    
    dt : int
        Time step used in integration in years.

    Returns
    -------
    2 x 2 panel comparing trajectories, total energy vs time, relative energy error vs time, and the virial ratio q over time for each method.

    '''
    import sys
    import os
    module_dir = os.path.abspath('../src') #directory to source

    sys.path.append(module_dir)
    from ODE import ODE
    import numpy as np
    import matplotlib.pyplot as plt

    #initial conditions to generate data
    ini_pos = np.array([[0, 0], [1, 0], [5.2, 0]])
    ini_vel = np.array([[0, 0], [0, 2*np.pi], [0, 0.88*np.pi]])
    masses = np.array([1, 3e-6, 9.55e-4])

    #ODE
    rk4 = ODE('rk4', 3, ini_pos, ini_vel, 2, duration, time_step, masses)
    leapfrog = ODE('leapfrog', 3, ini_pos, ini_vel, 2, duration, time_step, masses)
        
    #time
    time = np.linspace(0, duration, round(duration/time_step))
    
    #trajectories
    rk4_earth_x = rk4.positions[:,1,0]
    rk4_earth_y = rk4.positions[:,1,1]
    rk4_jupiter_x = rk4.positions[:,2,0]
    rk4_jupiter_y = rk4.positions[:,2,1]
    
    leapfrog_earth_x = leapfrog.positions[:,1,0]
    leapfrog_earth_y = leapfrog.positions[:,1,1]
    leapfrog_jupiter_x = leapfrog.positions[:,2,0]
    leapfrog_jupiter_y = leapfrog.positions[:,2,1]
    
    #relative energy error
    rk4_error = (rk4.E_tot - rk4.E_tot[0]) / np.abs(rk4.E_tot[0])
    leapfrog_error = (leapfrog.E_tot - leapfrog.E_tot[0]) / np.abs(leapfrog.E_tot[0])

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes[0, 0].plot(rk4_earth_x, rk4_earth_y, linestyle='-.', color='purple', label='RK4 earth')
    axes[0, 0].plot(leapfrog_earth_x, leapfrog_earth_y, linestyle='--', color='blue', label='Leapfrog earth')
    axes[0, 0].plot(rk4_jupiter_x, rk4_jupiter_y, linestyle='-.', color='red', label='RK4 jupiter')
    axes[0, 0].plot(leapfrog_jupiter_x, leapfrog_jupiter_y, linestyle='--', color='orange', label=' leapfrog jupiter')
    axes[0, 0].set_title('Trajectories')
    axes[0, 0].set_ylabel('X (AU)')
    axes[0, 0].set_xlabel('Y (AU)')
    axes[0, 0].grid(True)
    axes[0, 0].legend()     

    axes[0, 1].plot(time, rk4.E_tot, linestyle='-.', color='purple', label='RK4 total energy')
    axes[0, 1].plot(time, leapfrog.E_tot, linestyle='-.', color='red', label='Leapfrog total energy')
    axes[0, 1].plot(time, rk4.kE, linestyle='--', color='blue', label='RK4 kinetic energy')
    axes[0, 1].plot(time, leapfrog.kE, linestyle='--', color='green', label='Leapfrog kinetic energy')
    axes[0, 1].plot(time, rk4.U_g, linestyle=':', color='orange', label='RK4 gravitational potential energy')
    axes[0, 1].plot(time, leapfrog.U_g, linestyle=':', color='black', label='Leapfrog gravitational potential energy')
    axes[0, 1].set_title('Total Energy Change over Time')
    axes[0, 1].set_ylabel('Energy (M_solar AU^2/yr^2)')
    axes[0, 1].set_xlabel('Time (yrs)')
    axes[0, 1].grid(True)
    axes[0, 1].legend()    

    axes[1, 0].plot(time, rk4_error, color='purple', linestyle='-.',label='RK4')
    axes[1, 0].plot(time, leapfrog_error, color='red', linestyle=':', label='Leapfrog')
    #axes[1, 0].set_yscale('log')
    axes[1, 0].set_title('Relative Energy Error Change over Time')
    axes[1, 0].set_ylabel('Relative Energy Error')
    axes[1, 0].set_xlabel('Time (yrs)')
    axes[1, 0].grid(True)
    axes[1, 0].legend()    

    axes[1, 1].plot(time, rk4.virial, color='purple', linestyle='-.', label='RK4')
    axes[1, 1].plot(time, leapfrog.virial, color='red', linestyle=':', label='Leapfrog')
    axes[1, 1].set_title('Virial Change over Time')
    axes[1, 1].set_ylabel('Virial Ratio')
    axes[1, 1].set_xlabel('Time (yrs)')
    axes[1, 1].grid(True)
    axes[1, 1].legend()    

    fig.suptitle('Comparing Approximation Methods for 3 Body Earth Sun Problem')
    plt.tight_layout()
    plt.savefig('../outputs/figures/Comparing Approximation Methods for 3 Body Earth Sun Problem')
    plt.show()