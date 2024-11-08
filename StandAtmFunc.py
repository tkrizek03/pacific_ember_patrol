import math

def StdAtm(h):
    """
    Calculates standard atmosphere values.
    Altitude must be given in the troposphere or first isothermal region.

    Parameters
    ----------
    h : float
        Altitude in m.

    Returns
    -------
    h : float
        Geopotential altitude in m.
    t : float
        Temperature in K.
    d : float
        Density in kg/m**3.
    p : float
        Pressure in N/m**2.

    """
    
    # Sea Level Constants
    tSL = 288.16 # Temp; K
    dSL = 1.225 # Density; kg/m**3
    pSL = 1.01325e5 # Pressure; N/m**2
    
    # Other Constants
    g = 9.81 # Acceleration due to Gravity; m/s**2
    R = 287 # Gas Constant; J/kgK
    lapseRateT = -6.5e-3 # Lapse Rate in the Troposphere; K/m
    
    if h >= 0 and h < 11000:
        # Executes if geopotential altitude is within the troposphere.
        
        # Temperature
        t = tSL + (lapseRateT*h)
        
        # Pressure
        p = pSL * (t/tSL)**( (-g) / (lapseRateT*R) )
        
        # Density
        d = p / (R*t)
        
    elif h >= 11000 and h <= 25000:
        # Executes if geopotential altitude is within the isothermal layer.
        
        # Temperature
        t = tSL + (lapseRateT*11000)

        # Pressure
        pBase = pSL * (t/tSL)**( (-g) / (lapseRateT*R) )
        p = pBase * math.exp( ( (-g) / (R*t) ) * (h - 11000))
        
        # Density
        d = p / (R*t)
        
    return h, t, d, p