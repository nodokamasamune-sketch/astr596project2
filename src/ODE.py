import numpy as np

class ODE:
    '''
    Integrate gravitational interactions between N number of bodies via given approximation method.

    Parameters
    ----------
    method : str
        Integration method : euler, rk2, rk4, leapfrog
    
    N : float
        Number of bodies to calculate interactions for
    
    initial_positions : np.array
        Initial positions of every object in cartesian coordinates
        Shape N x dimensions
    
    initial_velocities : np.array
        Initial velocities of every object in cartesian coordinates
        Shape N x dimensions
    
    dimensions : float
        Number of dimensions to calculate interactions in.

    duration : int
        Number of years to simulate.
    
    dt : int
        Time step used in integration in years.
    
    masses : np.array
        Masses of every objection in solar masses, shape: N 

    R : int, opt
        Radius of cluster in AU.

    
    Returns
    -------
    Array of positions, velocities, and kinetic energy, potential energy, and virial ratio calculated at each step of the integration.
    positions : np.array
        Positions of every object in cartesian coordinates
        Shape nsteps+1 x N x dimensions 
     velocities : np.array
        Velocities of every object in cartesian coordinates
        Shape nsteps+1 x N x dimensions   
    kE : np.array
        Kinetic energy of system at every step of integration.
        Shape nsteps x dimensions
    U_g : np.array
        Gravitational potential energy of system at every step of integration.
        Shape nsteps x dimensions
    E_tot : np.array
        Total energy of system at every step of integration.
        Shape nsteps x dimensions
    virial : np.array
        Virial ratio of system at every step of integration.
        Shape nsteps x dimensions
    
    '''
    def __init__(self, method, N, initial_positions, initial_velocities, dimensions, duration, dt, masses=None, R=0):
        self.G = 4*np.pi**2
        if masses is not None:
            self.masses = masses

        self.N = N
        self.dimensions = dimensions
        self.R = R
            
        self.duration = duration
        self.dt = dt
        self.nsteps = round(duration/dt)
        
        self.positions = np.empty((self.nsteps+1, N, self.dimensions))
        self.velocities = np.empty((self.nsteps+1, N, self.dimensions))
        
        self.kE = np.empty(self.nsteps)
        self.U_g = np.empty(self.nsteps)
        self.E_tot = np.empty(self.nsteps)
        self.virial = np.empty(self.nsteps)
        
        
        self.positions[0,:] = initial_positions
        self.velocities[0,:] = initial_velocities

        self.e = 0.001 * self.R / self.N**(1/3)
        
        if method=='euler':
            for i in range(self.nsteps):
                new_positions, new_velocities, new_kE, new_U_g = self.euler(i)
                self.positions[i+1,:] = new_positions
                self.velocities[i+1,:] = new_velocities
                self.kE[i] = np.sum(new_kE)
                self.U_g[i] = np.sum(new_U_g) 
                if self.U_g[i] == 0:
                    print(self.U_g[i])
                self.E_tot[i] = np.sum(self.kE[i]) + np.sum(self.U_g[i])
                self.virial[i] = np.abs(2*self.kE[i] + self.U_g[i]) / np.abs(self.U_g[i])

        if method=='rk2':
            for i in range(self.nsteps):
                new_positions, new_velocities, new_kE, new_U_g = self.rk2(i)
                self.positions[i+1,:] = new_positions
                self.velocities[i+1,:] = new_velocities
                self.kE[i] = np.sum(new_kE)
                self.U_g[i] = np.sum(new_U_g) 
                if self.U_g[i] == 0:
                    print(self.U_g[i])
                self.E_tot[i] = np.sum(self.kE[i]) + np.sum(self.U_g[i])
                self.virial[i] = np.abs(2*self.kE[i] + self.U_g[i]) / np.abs(self.U_g[i])     

        if method=='rk4':
            for i in range(self.nsteps):
                new_positions, new_velocities, new_kE, new_U_g = self.rk4(i)
                self.positions[i+1,:] = new_positions
                self.velocities[i+1,:] = new_velocities
                self.kE[i] = np.sum(new_kE)
                self.U_g[i] = np.sum(new_U_g) 
                if self.U_g[i] == 0:
                    print(self.U_g[i])
                self.E_tot[i] = np.sum(self.kE[i]) + np.sum(self.U_g[i])
                self.virial[i] = np.abs(2*self.kE[i] + self.U_g[i]) / np.abs(self.U_g[i])

        if method=='leapfrog':
            for i in range(self.nsteps):
                new_positions, new_velocities, new_kE, new_U_g = self.leapfrog(i)
                self.positions[i+1,:] = new_positions
                self.velocities[i+1,:] = new_velocities
                self.kE[i] = np.sum(new_kE)
                self.U_g[i] = np.sum(new_U_g) 
                if self.U_g[i] == 0:
                    print(self.U_g[i])
                self.E_tot[i] = np.sum(self.kE[i]) + np.sum(self.U_g[i])
                self.virial[i] = np.abs(2*self.kE[i] + self.U_g[i]) / np.abs(self.U_g[i])


    def euler(self, i):
        positions = self.positions
        velocities = self.velocities
        N = self.N
        dt = self.dt
        self.i = i
        masses = self.masses
        
        current_positions = positions[self.i,:]
        current_velocities = velocities[self.i,:]
        
        new_positions = current_positions + current_velocities*dt

        accelerations = np.zeros_like(current_positions)
        new_U_g = np.zeros_like(current_positions)
        new_kE = np.zeros_like(current_positions)
        
        for i in range(N):
            new_kE = 0.5*masses[i]*np.linalg.norm(current_velocities[i])**2
            for j in range(N):
                if i != j:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij) 
                    accelerations[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)
                    #new_U_g = -self.G*masses[i]*masses[j]/r_mag
        
        for i in range(N):
            for j in range(N):
                if j > i:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij)
                    new_U_g = - (self.G * masses[i] * masses[j])/ (r_mag + self.e)

        new_velocities = current_velocities + accelerations*dt
        
        return new_positions, new_velocities, new_kE, new_U_g
            



    def rk2(self, i):
        positions = self.positions
        velocities = self.velocities
        N = self.N
        dt = self.dt
        self.i = i
        masses = self.masses
        
        current_positions = positions[self.i,:]
        current_velocities = velocities[self.i,:]
        
        accelerations = np.zeros_like(current_positions)
        new_U_g = np.zeros_like(current_positions)
        new_kE = np.zeros_like(current_positions)
        for i in range(N):
            new_kE = 0.5*masses[i]*np.linalg.norm(current_velocities[i])**2
            for j in range(N):
                if i != j:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij) 
                    accelerations[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)
                    #new_U_g = -self.G*masses[i]*masses[j]/r_mag

        for i in range(N):
            for j in range(N):
                if j > i:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij)
                    new_U_g = - (self.G * masses[i] * masses[j])/ (r_mag + self.e)

        vel_mid = current_velocities + 0.5 * accelerations * dt
        
        pos_mid = current_positions + 0.5 * current_velocities * dt

        accelerations_mid = np.zeros_like(current_positions)
        for i in range(N):
            for j in range(N):
                if i != j:
                    r_ij = pos_mid[j] - pos_mid[i]
                    r_mag = np.linalg.norm(r_ij)
                    accelerations_mid[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)
        

        new_positions = current_positions + vel_mid * dt
        new_velocities = current_velocities + accelerations_mid * dt
        return new_positions, new_velocities, new_kE, new_U_g
            
        
    def rk4(self, i):
        positions = self.positions
        velocities = self.velocities
        N = self.N
        dt = self.dt
        self.i = i
        masses = self.masses
        
        current_positions = positions[self.i,:]
        current_velocities = velocities[self.i,:]

        #step1
        p1 = current_positions
        v1 = current_velocities
        a1 = np.zeros_like(current_positions)
        new_U_g = np.zeros_like(current_positions)
        new_kE = np.zeros_like(current_positions)
        for i in range(N):
            new_kE = 0.5*masses[i]*np.linalg.norm(current_velocities[i])**2
            for j in range(N):
                if i != j:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij) 
                    a1[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)
                    #new_U_g = -self.G*masses[i]*masses[j]/r_mag

        for i in range(N):
            for j in range(N):
                if j > i:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij)
                    new_U_g = - (self.G * masses[i] * masses[j])/ (r_mag) + self.e

        #step2
        p2 = p1 + 0.5 * v1 * dt
        v2 = v1 + 0.5 * a1 * dt
        a2 = np.zeros_like(current_positions)
        for i in range(N):
            for j in range(N):
                if i != j:
                    r_ij = p2[j] - p2[i]
                    r_mag = np.linalg.norm(r_ij) 
                    a2[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)

        #step3
        p3 = p1 + 0.5 * v2 * dt
        v3 = v1 + 0.5 * a2 * dt
        a3 = np.zeros_like(current_positions)
        for i in range(N):
            for j in range(N):
                if i != j:
                    r_ij = p3[j]- p3[i]
                    r_mag = np.linalg.norm(r_ij) 
                    a3[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)        

        #step4
        p4 = p1 + v3 * dt
        v4 = v1 + a3 * dt
        a4 = np.zeros_like(current_positions)
        for i in range(N):
            for j in range(N):
                if i != j:
                    r_ij = p4[j]- p4[i]
                    r_mag = np.linalg.norm(r_ij) 
                    a4[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)           

        new_positions = current_positions + (dt / 6) * (v1 + 2*v2 + 2*v3 + v4)
        new_velocities = current_velocities + (dt / 6) * (a1 + 2*a2 + 2*a3 + a4)
        return new_positions, new_velocities, new_kE, new_U_g


    def leapfrog(self, i):
        positions = self.positions
        velocities = self.velocities
        N = self.N
        dt = self.dt
        self.i = i
        masses = self.masses
        
        current_positions = positions[self.i,:]
        current_velocities = velocities[self.i,:]

        accelerations = np.zeros_like(current_positions)
        new_U_g = np.zeros_like(current_positions)
        new_kE = np.zeros_like(current_positions)
        for i in range(N):
            new_kE = 0.5*masses[i]*np.linalg.norm(current_velocities[i])**2
            for j in range(N):
                if i != j:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij) 
                    accelerations[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)
                    #new_U_g = -self.G*masses[i]*masses[j]/r_mag

        for i in range(N):
            for j in range(N):
                if j > i:
                    r_ij = current_positions[j] - current_positions[i]
                    r_mag = np.linalg.norm(r_ij)
                    new_U_g = - (self.G * masses[i] * masses[j])/ (r_mag + self.e)

        v_mid = current_velocities + 0.5 * accelerations * dt
        new_positions = current_positions + v_mid * dt
        
        new_accelerations = np.zeros_like(current_positions)
        for i in range(N):
            for j in range(N):
                if i != j:
                    r_ij = new_positions[j] - new_positions[i]
                    r_mag = np.linalg.norm(r_ij) 
                    new_accelerations[i] += self.G *masses[j] * r_ij / (r_mag**2 + self.e**2)**(3/2)

        new_velocities = v_mid + 0.5 * new_accelerations * dt
        return new_positions, new_velocities, new_kE, new_U_g
        