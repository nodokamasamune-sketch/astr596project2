def cluster_animation(N, duration, time_step, num_frames=20):
    '''
    Creates animated gif of n bodies moving. Positions estimated using leapfrog integegration method.

    Parameters
    ----------
    N : int
        Number of bodies.
    
    Duration : float 
        Duration of simulation period, in years.

    time_step : float
        Size of step taken in integration, in years.

    num_frames : int, opt
        Number of frames in gif. Keep low.
    
    Returns
    -------
    An animated gif showing simulated stellar movement.
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
    
    output_dir = "animation_plots"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    def n_body_integrator(N, duration, time_step):
        masses = IMF(N)
        ini_pos, ini_vel, a = plummer(masses)
        le_fr = ODE('leapfrog', N, ini_pos, ini_vel, 3, duration, time_step, masses, R = 100)

        return le_fr

    lf = n_body_integrator(N, duration, time_step)

    #graph snapshots: 3d locations for t/t_total = 0, 0.25, 0.5, 0.75, 1
    
    nsteps = round(duration/time_step)

    num_frames = 20
    for i in range(num_frames):
        ti = lf.positions[i*round(1/num_frames * nsteps)] 

    for i in range(num_frames):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = ti[:,0]
        y = ti[:,1]
        z = ti[:,2]

        ax.scatter(x, y, z)
        ax.set_xlabel('x (AU)')
        ax.set_ylabel('y (AU)')
        ax.set_zlabel('z (AU)')
        #ax.set_xlim(np.mean(x) - 100, np.mean(x) + 100)
        #ax.set_ylim(np.mean(y) - 100, np.mean(y) + 100)
        #ax.set_zlim(np.mean(z) - 100, np.mean(z) + 100)
        plt.title('Cluster Evolution')
        filename = os.path.join(output_dir, f'plot_{i}.png')
        plt.savefig(filename)
        #plt.show()
        plt.close() #closes figure without displaying it. reduces clutter and frees memory.

    def gif(image_dir): #defining function to make gif
        from PIL import Image
        import os

        image_files = [f for f in os.listdir(image_dir) if f.endswith('.png')]

        images = [Image.open(os.path.join(image_dir, f)) for f in image_files]


        output_gif = '../outputs/figures/animation.gif'
        images[0].save(
            output_gif,
            save_all=True,
            append_images=images[1:],
            duration=100,  # ms between frames
            loop=0         # loop forever
        )

        print(f"GIF saved as: {output_gif}")
    
    gif(output_dir) #calling function to make gif

def cluster_artificial(N, duration, time_step):
     '''
   Creates animated gif of artificial n bodies moving. Positions estimated using leapfrog integegration method. Done to produce 'more exciting'        animated gif.

    Parameters
    ----------
    N : int
        Number of bodies.
    
    Duration : float 
        Duration of simulation period, in years.

    time_step : float
        Size of step taken in integration, in years.

    num_frames : int, opt
        Number of frames in gif. Keep low.
    
    Returns
    -------
    An animated gif showing simulated stellar movement.   
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
    
    output_dir = "animation_plots_artificial"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
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

    ini_pos, ini_vel, masses = star_generator_random(N) #N = 10
        
    le_fr = ODE('leapfrog', N, ini_pos, ini_vel, 3, duration, time_step, masses, R = 10)

    #graph snapshots: 3d locations for t/t_total = 0, 0.25, 0.5, 0.75, 1
    
    nsteps = round(duration/time_step)

    num_frames = 40
    t = {}
    for i in range(num_frames):
        t[i] = le_fr.positions[i*round(1/num_frames * nsteps)] 

    for i in range(num_frames):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = t[i][:,0]
        y = t[i][:,1]
        z = t[i][:,2]

        ax.scatter(x, y, z)
        ax.set_xlabel('x (AU)')
        ax.set_ylabel('y (AU)')
        ax.set_zlabel('z (AU)')
        ax.set_xlim(-4000, 4000)
        ax.set_ylim(-4000, 4000)
        ax.set_zlim(-4000, 4000)
        #plt.title('Cluster at t/duration = 0')
        #plt.savefig('animation/randcluster',i,'.png')
        filename = os.path.join(output_dir, f'artificial_plot_{i}.png')
        plt.savefig(filename)
        #plt.show()
        plt.close() 
    def gif(image_dir): #defining function to make gif
        from PIL import Image
        import os

        image_files = [f for f in os.listdir(image_dir) if f.endswith('.png')]

        images = [Image.open(os.path.join(image_dir, f)) for f in image_files]


        output_gif = '../outputs/figures/animation_artificial.gif'
        images[0].save(
            output_gif,
            save_all=True,
            append_images=images[1:],
            duration=500,  # ms between frames
            loop=0         # loop forever
        )

        print(f"GIF saved as: {output_gif}")
    
    gif(output_dir) #calling function to make gif
    