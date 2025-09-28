def earth_sun_plotter():
    '''
    Plotting 2 body earth sun interaction. Produces :
        3 panel figure showing energy components vs. time, relative energy error vs. time, and virial ratio q over time for each method : euler,  RK2, RK4, and leapfrog. 
    2 x 2 panel comparing trajectories, total energy vs time, relative energy error vs time, and the virial ratio q over time for each method.

    Parameters
    ----------
    None: all parameters preset (such as initial positions of objects, velocities)
    
    Returns
    -------
    3 panel figure showing energy components vs. time, relative energy error vs. time, and virial ratio q over time for each method : euler,  RK2, RK4, and leapfrog. 
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
    ini_pos = np.array([[0, 0], [1, 0]])
    ini_vel = np.array([[0, 0], [0, 2*np.pi]])
    masses = np.array([1, 3e-6])

    #ODE- integrate from initial state using all four methods
    euler = ODE('euler', 2, ini_pos, ini_vel, 2, 10, 0.01, masses)
    rk2 = ODE('rk2', 2, ini_pos, ini_vel, 2, 10, 0.01, masses)
    rk4 = ODE('rk4', 2, ini_pos, ini_vel, 2, 10, 0.01, masses)
    leapfrog = ODE('leapfrog', 2, ini_pos, ini_vel, 2, 10, 0.01, masses)

    #data to plot

    #time 
    time = np.linspace(0, 10, 1000)
    
    #trajectories- split into x and y components for each method

    #euler
    euler_x = euler.positions[:,1:,0]
    euler_y = euler.positions[:,1,1]

    #rk2
    rk2_x = rk2.positions[:,1:,0]
    rk2_y = rk2.positions[:,1,1]

    #rk4
    rk4_x = rk4.positions[:,1:,0]
    rk4_y = rk4.positions[:,1,1]

    #leapfrog
    leapfrog_x = leapfrog.positions[:,1:,0]
    leapfrog_y = leapfrog.positions[:,1,1]
    
    
    #relative energy error = (E_tot - E_tot_i) / np.abs(E_tot_i)
    euler_error = (euler.E_tot-euler.E_tot[1]) / np.abs(euler.E_tot[1]) #euler
    rk2_error = (rk2.E_tot - rk2.E_tot[1]) / np.abs(rk2.E_tot[1]) #rk2
    rk4_error = (rk4.E_tot - rk4.E_tot[1]) / np.abs(rk4.E_tot[1]) #rk4
    leapfrog_error = (leapfrog.E_tot - leapfrog.E_tot[1]) / np.abs(leapfrog.E_tot[1]) #leapfrog
    
    
    #plotting
    
    #3 panels per method

    #euler
    #define figure, axes, and figsize
    fig, axes = plt.subplots(3, 1, figsize=(5,7))

    #first plot: kinetic energy, total energy,and potential energy vs. time
    axes[0].plot(time, euler.kE, color='blue', label='Kinetic Energy')
    axes[0].plot(time, euler.E_tot, color='orange', label='Total Energy')
    axes[0].plot(time, euler.U_g, color='purple', label='Gravitational Potential Energy')
    axes[0].set_title('Energy Change over Time')
    axes[0].set_ylabel('Energy (M_s * AU^2/ yr^2)')
    axes[0].set_xlabel('Time (yrs)')
    axes[0].grid(True)
    axes[0].legend()

    #second plot: relative energy error (log) vs. time
    axes[1].plot(time, euler_error, color='blue', label='Relative Energy Error')
    axes[1].set_yscale('log')
    axes[1].set_title('Relative Energy Error Change over Time')
    axes[1].set_ylabel('Relative Energy Error')
    axes[1].set_xlabel('Time (yrs)')
    axes[1].grid(True)
    axes[1].legend()    

    #third plot: virial ratio = |2*kE + U_g| / U_g vs. time
    axes[2].plot(time, euler.virial, color='blue', label='Virial')
    axes[2].set_title('Virial Change over Time')
    axes[2].set_ylabel('Virial Ratio')
    axes[2].set_xlabel('Time (yrs)')
    axes[2].grid(True)
    axes[2].legend()    

    fig.suptitle('Euler Method Approximation')
    plt.tight_layout()
    plt.savefig('../outputs/figures/2 Body Euler Approximation.png')
    #plt.close() #make figure and save as .png without displaying it
    plt.show()

    #rk2
    #define figure, axes, and figsize
    fig, axes = plt.subplots(3, 1, figsize=(5,7))

    #first plot: kinetic energy, total energy,and potential energy vs. time
    axes[0].plot(time, rk2.kE, color='blue', label='Kinetic Energy')
    axes[0].plot(time, rk2.E_tot, color='orange', label='Total Energy')
    axes[0].plot(time, rk2.U_g, color='purple', label='Gravitational Potential Energy')
    axes[0].set_title('Energy Change over Time')
    axes[0].set_ylabel('Energy (M_s * AU^2/ yr^2)')
    axes[0].set_xlabel('Time (yrs)')
    axes[0].grid(True)
    axes[0].legend()

    #second plot: relative energy error (log) vs. time
    axes[1].plot(time, rk2_error, color='blue', label='Relative Energy Error')
    axes[1].set_yscale('log')
    axes[1].set_title('Relative Energy Error Change over Time')
    axes[1].set_ylabel('Relative Energy Error')
    axes[1].set_xlabel('Time (yrs)')
    axes[1].grid(True)
    axes[1].legend()    

    #third plot: virial ratio = |2*kE + U_g| / U_g vs. time
    axes[2].plot(time, rk2.virial, color='blue', label='Virial')
    axes[2].set_title('Virial Change over Time')
    axes[2].set_ylabel('Virial Ratio')
    axes[2].set_xlabel('Time (yrs)')
    axes[2].grid(True)
    axes[2].legend()    

    fig.suptitle('RK2 Method Approximation')
    plt.tight_layout()
    plt.savefig('../outputs/figures/2 Body RK2 Approximation.png')
    #plt.close() #make figure and save as .png without displaying it
    plt.show()

    #RK4
    fig, axes = plt.subplots(3, 1, figsize=(5,7))
    
    #first plot: kinetic energy, total energy,and potential energy vs. time
    axes[0].plot(time, rk4.kE, color='blue', label='Kinetic Energy')
    axes[0].plot(time, rk4.E_tot, color='orange', label='Total Energy')
    axes[0].plot(time, rk4.U_g, color='purple', label='Gravitational Potential Energy')
    axes[0].set_title('Energy Change over Time')
    axes[0].set_ylabel('Energy (M_s * AU^2/ yr^2)')
    axes[0].set_xlabel('Time (yrs)')
    axes[0].grid(True)
    axes[0].legend()

    #second plot: relative energy error (log) vs. time
    axes[1].plot(time, rk4_error, color='blue', label='Relative Energy Error')
    axes[1].set_yscale('log')
    axes[1].set_title('Relative Energy Error Change over Time')
    axes[1].set_ylabel('Relative Energy Error')
    axes[1].set_xlabel('Time (yrs)')
    axes[1].grid(True)
    axes[1].legend()    

    #third plot: virial ratio = |2*kE + U_g| / U_g vs. time
    axes[2].plot(time, rk4.virial, color='blue', label='Virial')
    axes[2].set_title('Virial Change over Time')
    axes[2].set_ylabel('Virial')
    axes[2].set_xlabel('Time (yrs)')
    axes[2].grid(True)
    axes[2].legend()    

    fig.suptitle('RK4 Method Approximation')
    plt.tight_layout()
    plt.savefig('../outputs/figures/2 Body RK4 Approximation.png')
    #plt.close() #make figure and save as .png without displaying it
    plt.show()

    #leapfrog
    fig, axes = plt.subplots(3, 1, figsize=(5,7))

    #first plot: kinetic energy, total energy,and potential energy vs. time
    axes[0].plot(time, leapfrog.kE, color='blue', label='Kinetic Energy')
    axes[0].plot(time, leapfrog.E_tot, color='orange', label='Total Energy')
    axes[0].plot(time, leapfrog.U_g, color='purple', label='Gravitational Potential Energy')
    axes[0].set_title('Energy Change over Time')
    axes[0].set_ylabel(('Energy (M_s * AU^2/ yr^2)'))
    axes[0].set_xlabel('Time (yrs)')
    axes[0].grid(True)
    axes[0].legend()

    #second plot: relative energy error (log) vs. time
    axes[1].plot(time, leapfrog_error, color='blue', label='Relative Energy Error')
    axes[1].set_yscale('log')
    axes[1].set_title('Relative Energy Error Change over Time')
    axes[1].set_ylabel('Relative Energy Error')
    axes[1].set_xlabel('Time (yrs)')
    axes[1].grid(True)
    axes[1].legend()    

    #third plot: virial ratio = |2*kE + U_g| / U_g vs. time
    axes[2].plot(time, leapfrog.virial, color='blue', label='Virial')
    axes[2].set_title('Virial Change over Time')
    axes[2].set_ylabel('Virial')
    axes[2].set_xlabel('Time (yrs)')
    axes[2].grid(True)
    axes[2].legend()    

    fig.suptitle('Leapfrog Method Approximation')
    plt.tight_layout()
    plt.savefig('../outputs/figures/2 Body Leapfrog Approximation.png')
    #plt.close() #make figure and save as .png without displaying it
    plt.show()


    #2x2 comparison figure for all four methods showing:
    #orbital tragectories - different line styles
    #total energy vs time - different colors
    #rel energy error - log scale
    #virial ratio
    
    #trajectories
    euler_x = euler.positions[:,1:,0]
    euler_y = euler.positions[:,1,1]

    rk2_x = rk2.positions[:,1:,0]
    rk2_y = rk2.positions[:,1,1]

    rk4_x = rk4.positions[:,1:,0]
    rk4_y = rk4.positions[:,1,1]
    
    leapfrog_x = leapfrog.positions[:,1:,0]
    leapfrog_y = leapfrog.positions[:,1,1]
        
    #defining 2 x 2 subplot + figsize
    fig, axes = plt.subplots(2, 2, figsize=(6, 6))

    #orbital tragectories - different line styles
    axes[0, 0].plot(euler_x, euler_y, linestyle='-',color='blue', label='Euler')
    axes[0, 0].plot(rk2_x, rk2_y, linestyle='--', color='orange', label='RK2')
    axes[0, 0].plot(rk4_x, rk4_y, linestyle='-.', color='purple', label='RK4')
    axes[0, 0].plot(leapfrog_x, leapfrog_y, linestyle=':', color='grey', label='Leapfrog')
    axes[0, 0].set_title('Trajectories')
    axes[0, 0].set_ylabel('X (AU)')
    axes[0, 0].set_xlabel('Y (AU)')
    axes[0, 0].grid(True)
    axes[0, 0].legend()     

    #total energy vs time - different colors
    axes[0, 1].plot(time, euler.E_tot, linestyle='--',color='blue', label='Euler')
    axes[0, 1].plot(time, rk2.E_tot, linestyle=':', color='orange', label='RK2')
    axes[0, 1].plot(time, rk4.E_tot, linestyle='-.', color='purple', label='RK4')
    axes[0, 1].plot(time, leapfrog.E_tot, linestyle=':', color='red', label='Leapfrog')
    axes[0, 1].set_title('Total Energy Change over Time')
    axes[0, 1].set_ylabel(('Energy (M_s * AU^2/ yr^2)'))
    axes[0, 1].set_xlabel('Time (yrs)')
    axes[0, 1].grid(True)
    axes[0, 1].legend()    

    #rel energy error - log scale
    axes[1, 0].plot(time, euler_error, linestyle='--', color='blue', label='Euler')
    axes[1, 0].plot(time, rk2_error, linestyle=':', color='orange', label='RK2')
    axes[1, 0].plot(time, rk4_error, linestyle='-.',  color='purple', label='RK4')
    axes[1, 0].plot(time, leapfrog_error, linestyle=':', color='red', label='Leapfrog')
    axes[1, 0].set_yscale('log')
    axes[1, 0].set_title('Relative Energy Error Change over Time')
    axes[1, 0].set_ylabel('Relative Energy Error')
    axes[1, 0].set_xlabel('Time (yrs)')
    axes[1, 0].grid(True)
    axes[1, 0].legend()    

    #virial ratio
    axes[1, 1].plot(time, euler.virial, linestyle='--', color='blue', label='Euler')
    axes[1, 1].plot(time, rk2.virial, linestyle=':', color='orange', label='RK2')
    axes[1, 1].plot(time, rk4.virial, linestyle='-.', color='purple', label='RK4')
    axes[1, 1].plot(time, leapfrog.virial, linestyle=':', color='red', label='Leapfrog')
    axes[1, 1].set_title('Virial Change over Time')
    axes[1, 1].set_ylabel('Virial')
    axes[1, 1].set_xlabel('Time (yrs)')
    axes[1, 1].grid(True)
    axes[1, 1].legend()    

    fig.suptitle('Comparing Approximation Methods for 2 Body Earth Sun Problem')
    plt.tight_layout()
    plt.savefig('../outputs/figures/Comparing Approximation Methods for 2 Body Earth Sun Problem')
    #plt.close() #make figure and save as .png without displaying it
    plt.show()