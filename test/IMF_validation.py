#IMF_validation

#checking P1/P2 proportions match outputed masses

def proportion_check(N):
    '''
    Checking that masses produced by IMF follows distribution set by IMF. Prints:
        - p1, p2, and mass ratio below m_b.
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
    from IMF_plummer import IMF

    masses = IMF(N)

    over = masses > 0.5

    #print(len(over)/N)
    #print(np.sum(over))
    print(1 - np.sum(over)/N)