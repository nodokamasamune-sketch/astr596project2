def cluster_evolution_plotter(N, duration, time_step):
    '''
    Plotting N body interaction. Produces :
        5 panel figure showing various stages of cluster evolution. Positions calculated using leapfrog approximation.

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
    5 panel figure showing various stages of cluster evolution.
    '''
    import sys
    import os
    module_dir = os.path.abspath('../src') #directory to source

    sys.path.append(module_dir)
    import numpy as np
    import matplotlib.pyplot as plt
    from ODE import ODE
    from IMF_plummer import IMF, plummer
    from mpl_toolkits import mplot3d
    
    def n_body_integrator(N, duration, time_step):
        masses = IMF(N)
        ini_pos, ini_vel, a = plummer(masses)
        lf = ODE('leapfrog', N, ini_pos, ini_vel, 3, duration, time_step, masses, R = 100)

        return lf

    lf = n_body_integrator(N, duration, time_step)

    #graph snapshots: 3d locations for t/t_total = 0, 0.25, 0.5, 0.75, 1

    nsteps = round(duration/time_step)

    t0 = lf.positions[0]
    t1 = lf.positions[round(0.25*nsteps)]
    t2 = lf.positions[round(0.5*nsteps)]
    t3 = lf.positions[round(0.75*nsteps)]
    t4 = lf.positions[-1]


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = t0[:,0]
    y = t0[:,1]
    z = t0[:,2]

    ax.scatter(x, y, z)
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_zlabel('z (AU)')
    plt.title('Cluster at t/duration = 0')
    plt.savefig('../outputs/figures/Cluster at t = 0.png')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = t1[:,0]
    y = t1[:,1]
    z = t1[:,2]

    ax.scatter(x, y, z)
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_zlabel('z (AU)')
    plt.title('Cluster at t/duration = 0.25')
    plt.savefig('../outputs/figures/Cluster at t = quarter way through.png')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = t2[:,0]
    y = t2[:,1]
    z = t2[:,2]

    ax.scatter(x, y, z)
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_zlabel('z (AU)')
    plt.title('Cluster at t/duration = 0.5')
    plt.savefig('../outputs/figures/Cluster at t = halfway through.png')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = t3[:,0]
    y = t3[:,1]
    z = t3[:,2]

    ax.scatter(x, y, z)
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_zlabel('z (AU)')
    plt.title('Cluster at t/duration = 0.75')
    plt.savefig('../outputs/figures/Cluster at t = three quarters way through.png')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = t4[:,0]
    y = t4[:,1]
    z = t4[:,2]

    ax.scatter(x, y, z)
    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_zlabel('z (AU)')
    plt.savefig('../outputs/figures/Cluster at final point.png')
    plt.title('Cluster at t/duration = 1')
    plt.show()