"""

Help keyword: Numpy array slicing
    
    Syntax [iZ,:] means, we are accessing the whole row iZ
    Syntax [:,iX] means, we are accessing the whole column iX
    Syntax [:idZ,iX] means, we're accessing elements 0..idZ-1 in column iX
    Syntax [idZ:,iX] means, we're accessing elements idZ:end in column iX
    
    Compare e.g.
      https://www.programiz.com/python-programming/numpy/array-slicing
    or
      https://www.w3schools.com/python/numpy/numpy_array_slicing.asp


"""

import numpy as np
from mpmath import coth



#%% Function definitions
def calcSurfElevAiry(ampli, theta):
    '''
    Calculate surface elevation at point(s) xEval.

    Parameters
    ----------
    ampli : float
        amplitude
    theta : float
        wave phase

    Returns
    -------
    eta : float
        surface elevation

    '''
    eta = ampli*np.cos(theta)
    
    return eta
    
    
def calcVelAiry(ampli, kNum, om, phi, depth, xEval, zEval, tEval):
    '''
    Calculate velocity at point xEval, zEval.
    
    Parameters
    ----------
    ampli : float
        amplitude
    kNum : float
        wave-number k
    om : float
        angular frequency
    depth : float
        depth
    xEval : float
        x-position
    zEval : float
        z-position

    Returns
    -------
    uX : float
        velocity field x-component
    uZ : float
        velocity field x-component

    '''
    # Needed to catch numerically ill conditioned deep water conditions
    import warnings
    warnings.filterwarnings("error")
    
    
    # Phase angle
    theta  = kNum*xEval - om*tEval + phi
    
    
    # General solution
    try:
        # Constant factors
        ampliC = om*ampli/np.sinh(kNum*depth)
    
        # Depth dependent parameters
        cosh_kzd = np.cosh(kNum*(zEval+depth))
        sinh_kzd = np.sinh(kNum*(zEval+depth))
        
        # Velocity field
        uX = ampliC*(cosh_kzd)*np.cos(theta)
        uZ = ampliC*(sinh_kzd)*np.sin(theta)
    
    
    # Deep water solution
    except RuntimeWarning:
        uX = om*ampli*np.exp(kNum*zEval)*np.cos(theta)
        uZ = om*ampli*np.exp(kNum*zEval)*np.sin(theta)
    
    
    # Reset to normal warning behaviour
    warnings.resetwarnings()
    
    return uX, uZ



def calcVelFieldAiry(ampli, kNum, om, phi, depth, xEval, zEval, tEval):
    '''
    Calculate velocity field and surface elevation on coordinate grid.
    Slightly improved version making use of numpy vectorisation as compared
    to iterating over list of single points.
    
    Parameters
    ----------
    ampli : float
        amplitude
    kNum : float
        wave-number k
    om : float
        angular frequency
    depth : float
        depth
    xEval : 1d numpy array
        x-positions
    zEval : 1d numpy array
        z-positions

    Returns
    -------
    eta : 1d numpy array
        surface elevation
    uX : 2d numpy array
        velocity field x-component
    uZ : 2d numpy array
        velocity field x-component

    '''
    # Needed to catch numerically ill conditioned deep water conditions
    import warnings
    warnings.filterwarnings("error")
    
    
    # Allocate output arrays
    nPtsX = xEval.shape[0]
    nPtsZ = zEval.shape[0]

    uX = np.zeros((nPtsZ,nPtsX))
    uZ = np.zeros((nPtsZ,nPtsX))
    
    
    # Phase angles
    theta  = kNum*xEval - om*tEval + phi
    stheta = np.sin(theta)
    ctheta = np.cos(theta)
    
    # Surface elevation
    eta = calcSurfElevAiry(ampli, theta)
    
    
    # General solution
    try:
        # Constant factors
        ampliC = om*ampli/np.sinh(kNum*depth)
        
        # Depth dependent parameters
        cosh_kzd = np.cosh(kNum*(zEval+depth))
        sinh_kzd = np.sinh(kNum*(zEval+depth))
        
        # Velocity field
        for iZ in range(nPtsZ):
            uX[iZ,:] = ampliC*cosh_kzd[iZ]*ctheta
            uZ[iZ,:] = ampliC*sinh_kzd[iZ]*stheta
    
    
    # Deep water solution
    except RuntimeWarning:
        eTokz = np.exp(kNum*zEval)
        
        for iZ in range(nPtsZ):
            uX[iZ,:] = om*ampli*eTokz[iZ]*ctheta
            uZ[iZ,:] = om*ampli*eTokz[iZ]*stheta
        
    
    # Reset to normal warning behaviour
    warnings.resetwarnings()
    
    return eta, uX, uZ



def clipValuesAboveSurface(eta, uX, uZ, zEval):
    '''
    Clips values to numpy.nan ("not a number") above surface.

    Parameters
    ----------
    eta : 1d numpy array
        surface elevation
    uX : 2d numpy array
        velocity field x-component
    uZ : 2d numpy array
        velocity field x-component
    zEval : 1d numpy array
        z-positions

    Returns
    -------
    uX : 2d numpy array
        velocity field x-component
    uZ : 2d numpy array
        velocity field x-component

    '''
    nPtsX = uX.shape[1]
    dZ    = (zEval[1]-zEval[0])
    
    # We are always counting down from the highest z-values to the surface
    # Case 1: zEval stored in descending order
    if dZ < 0:
        for iX in range(nPtsX):
            idZ = int((eta[iX] - zEval[0])/dZ)
            
            uX[:idZ,iX] = np.nan
            uZ[:idZ,iX] = np.nan
        
    # Case 2: zEval stored in ascending order
    else:
        for iX in range(nPtsX):
            idZ = int((eta[iX] - zEval[0])/dZ)
            
            uX[idZ:,iX] = np.nan
            uZ[idZ:,iX] = np.nan
    
        
                
    return uX, uZ



#%% Spectra and wave lengths
def calcSpectrumBretschneider( fMax, deltaF, HSig, TPeak, depth, gAcc, nIter=10 ):
    frequs  = np.linspace(deltaF,fMax,round(fMax/deltaF))
    oms     = 2*np.pi*frequs
    deltaOm = 2.0*np.pi*deltaF
    
    periods = 1.0/frequs
    nWaves  = len(periods)

    spectrum = np.zeros(nWaves)
    
    omPeak = 2*np.pi/TPeak
    
    spectrum = 5.0/16.0*(HSig**2)*(omPeak**4)*(oms**-5)* \
               np.exp(-5.0/4.0*(omPeak/oms)**4)        

    amps = np.sqrt(2*spectrum*deltaOm)
    
    kNumbers = np.zeros(nWaves)
    wLengths = np.zeros(nWaves)

    for iW in range(nWaves):
        kNumbers[iW] = calcWavenumberNewt(depth, oms[iW], gAcc, nIter)
    
    wLengths = 2*np.pi/kNumbers
             
    return nWaves, periods, oms, spectrum, amps, wLengths, kNumbers


def calcWavenumberNewt(h,om,g,nIter):
    k0  = om**2/g
    k0h = k0*h
    
    # Initial guess
    if k0h >= 1.0:
        kh = k0*h
    else:
        kh = (k0*h)**0.5
    
    for i in range(nIter):
        cothkh = coth(kh)
        kh = kh - (kh-k0h*cothkh)/(1.0+k0h*(cothkh**2-1.0))
    
    return (kh/h)


def calcWavePeriod(h,L,g):
    k    = 2*np.pi/L
    omSq = g*k*np.tanh(k*h)
    
    return 2*np.pi/np.sqrt(omSq)

