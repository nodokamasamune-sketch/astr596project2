def IMF(N):
    '''
    Generate masses using Kroupa Initial Mass Function.

    Parameters
    ----------
    N : float
        Number of bodies to calculate masses for

    Returns
    -------
    masses : np.array
        Masses, shape N
    '''
    import numpy as np
    
    #constants
    alpha1 = 1.3
    alpha2 = 2.3
    m_min = 0.08
    m_b = 0.5
    m_max = 150


    #precompute probabilities p1 & p2
    b_1 = m_b**(1- alpha1)
    a_1 = m_min**(1-alpha1)
    b_2 = m_max**(1-alpha2)
    a_2 = m_b**(1-alpha2)
    
    #I1 = (m_b**(1 - alpha1) - m_min**(1 - alpha1))/ (1 - alpha1)
    #I2 = (m_max**(1 - alpha2) - m_b**(1 - alpha2))/ (1 - alpha2)

    I1 = (b_1 - a_1) / (1 - alpha1)
    I2 = (b_2 - a_2) / (1 - alpha2)
    
    A1 = N / (I1 + I2*(m_b**(alpha2 - alpha1)))
    A2 = A1*(m_b)**(alpha2 - alpha1)

    p1 = (A1*I1)/N
    p2 = (A2*I2)/N

    #print(p1)
    print(p2 + p1)

    #generate probability
    #U = np.random.uniform(10**-12, 1, N)
    
    masses = np.zeros(N)
    
    for i in range(N):
        u = np.random.uniform(10**-12, 1)
        s = np.random.uniform(0, 1)
        
        if s < p1:
            m = (a_1 + u*(b_1 - a_1))**(1/(1-alpha1))
            
        else:
            m = (a_2 + u*(b_2 - a_2))**(1/(1-alpha2))

        masses[i] = m

    return masses

def plummer(masses):
    '''
    Generate initial positions and velocities given masses using Plummer sphere radial density model. Positions and velocities are corrected so that their mass weighted average is zero.

    Parameters
    ----------
    masses : np.array
        Masses, shape N

    Returns
    -------
    initial_positions : np.array 
        Initial position for each star in cartesian coordinates, shape N x 3
    initial_velocities : np.array
        Initial velocity for each star in cartesian coordinates, shape N x 3
    a : int
        Scaling factor used in calculations.
    '''
    import numpy as np
    N = len(masses)

    if N < 100:
        a = 100
    else:
        a = 1000
    
    coordinates = np.zeros([N, 3])
    
    for i in range(N):
        u = np.random.uniform(0, 1)
        r = a * (u**(-2/3) - 1)**(-1/2)

        u_phi = np.random.uniform(0, 1)
        u_theta = np.random.uniform(0, 1)
        phi = 2*np.pi*u_phi
        cos_theta = 1 - 2*u_theta

        x = r * np.sin(np.arccos(cos_theta)) * np.cos(u_phi)
        y = r * np.sin(np.arccos(cos_theta)) * np.sin(u_phi)
        z = r * cos_theta

        coords = np.array([x, y, z])
        #coordinates[i:0] = x
        #coordinates[i:1] = y
        #coordinates[i:2] = z

        coordinates[i,:] = coords

    x_avg = np.sum(coordinates[:,0] * masses) / np.sum(masses)
    y_avg = np.sum(coordinates[:,1] * masses) / np.sum(masses)
    z_avg = np.sum(coordinates[:,2] * masses) / np.sum(masses)

    center_correction = np.array([x_avg, y_avg, z_avg])
    centered_coords = coordinates - center_correction
    
    '''
    x_corr = np.sum(centered_coords[:,0,0] * masses)
    y_corr = np.sum(centered_coords[:,0,1] * masses)
    z_corr = np.sum(centered_coords[:,0,2] * masses)

    print(x_corr, y_corr, z_corr)
    '''

    velocities = np.zeros([N, 3])
    for i in range(N):
        dist = ((4*np.pi**2)/(6*a) *(1 + (np.linalg.norm(centered_coords[i]))**2 / a**2)**(-1/2)) ** (1/2)
        #print(dist)
        
        vel_dis = np.random.uniform(0, dist, 100)
        v_x = np.random.choice(vel_dis)
        v_y = np.random.choice(vel_dis)
        v_z = np.random.choice(vel_dis)
        #print(v_x, v_y, v_z)
        v_avg = ((v_x**2 + v_y**2 + v_z**2)/3 )**(1/2)
        dist_check = 3**(1/2)*dist
        #print(v_avg, dist_check)

        
        vel = np.array([v_x, v_y, v_z])
        velocities[i,:] = vel

    vx_avg = np.sum(velocities[:,0] * masses) / np.sum(masses)
    vy_avg = np.sum(velocities[:,1] * masses) / np.sum(masses)
    vz_avg = np.sum(velocities[:,2] * masses) / np.sum(masses)  

    center_velocities = np.array([vx_avg, vy_avg, vz_avg])
    centered_velocities = velocities - center_velocities
    
    vx_corr = np.sum(centered_velocities[:,0] * masses)
    vy_corr = np.sum(centered_velocities[:,1] * masses)
    vz_corr = np.sum(centered_velocities[:,2] * masses)

    #print(vx_corr, vy_corr, vz_corr)

    return centered_coords, centered_velocities, a

    