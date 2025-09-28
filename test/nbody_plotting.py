def nbody_plotting(N, duration, time_step):
    '''
    Produces plots of energy vs. time, relative energy error, and the virial ratio vs time 
    for an N body problem integrated using the leapfrog method. 

    Parameters
    ----------
    N : int
        Number of bodies

    duration : float
        Duration of simulation, in years

    time_step : float
        Size of step taken in integration process, in years.

    Returns
    -------
    3 panel figure with plots of kinetic energy, gravitational potential energy, and total energy vs. time,
    the relative energy error vs. time, and the virial ratio vs. time.
    '''
    import sys
    import os
    module_dir = os.path.abspath('../src') #directory to source
    sys.path.append(module_dir)
    import numpy as np
    import matplotlib.pyplot as plt
    from ODE import ODE
    from IMF_plummer import IMF, plummer
    
    masses = IMF(N)
    ini_pos, ini_vel, a = plummer(masses)
    lf = ODE('leapfrog', N, ini_pos, ini_vel, 3, duration, time_step, masses, R = 100)

    #calculate relative error plot
    leapfrog_error = (lf.E_tot - lf.E_tot[1]) / np.abs(lf.E_tot[1])

    #time in years
    time = np.linspace(0, duration, round(duration/time_step))

    fig, axes = plt.subplots(3, 1, figsize=(8, 10))

    #graph of kE, W, and E_tot on same axes
    axes[0].plot(time, lf.kE, color='red', label='Kinetic Energy')
    axes[0].plot(time, lf.U_g, color='green', label='Gravitational Potential Energy')
    axes[0].plot(time, lf.E_tot, color='blue', label='Total Energy')
    axes[0].set_ylabel('Energy (M_solar AU^2/yr^2)')
    axes[0].set_xlabel('Time (years)')
    axes[0].legend()

    #relative energy error over time
    axes[1].plot(time, leapfrog_error, color='red', label='Relative Error')
    axes[1].set_ylabel('Relative Energy Error')
    #axes[1].set_yscale('log')
    axes[1].set_xlabel('Time (years)')
    axes[1].legend()

    #virial ratio vs. time.
    axes[2].plot(time, lf.virial, color= 'purple', label='Virial')
    axes[2].set_ylabel('Virial Ratio')
    axes[2].set_xlabel('Time (years)')
    axes[2].legend()    

    fig.suptitle('N body simulation using the leapfrog integration method')
    plt.savefig('../outputs/figures/N body simulation using the leapfrog integration method')
    plt.tight_layout()
    plt.show()