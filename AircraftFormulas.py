import math

"""
Pacific Ember Patrol

Class for aircraft and flight condition advanced calcualtions.
"""

#############################

class AircraftFormulas:

    ###
    
    def coefficientDragComponent(cfc, ff, qc, sWet, sWing):
        """
        Calculates the component drag coefficient

        Parameters
        ----------
        cfc : float
            Component flat plate skin friction coefficient.
        ff : float
            Component form factor.
        qc : float
            Component interference factor.
        sWet : float
            Wetted area of the component.
        sWing : float
            Planform area of the wing

        Returns
        -------
        cD0component : float
            Component drag coefficient (unitless).

        """
        cD0component = ( cfc * ff * qc * sWet ) / (sWing)
            
        return cD0component
    
    ###
    

    def dynamicViscosity(t):
        """
        Calculates the dynamic viscosity.

        Parameters
        ----------
        t : float
            Temperature in kelvin.

        Returns
        -------
        mu : float
            mu value in (Ns)/m**2.

        """
        
        coef1 = 1.716e-5
        coef2 = (t/273.15)**1.5
        coef3 = ((273.15 + 113) / (t + 113))
        mu = coef1 * coef2 * coef3
        
        return mu
            
    ###
    
    def mach(temp, velocity):
        """
        Calculates mach number

        Parameters
        ----------
        temp : float
            Temperature in kelvin.
        velocity : float
            Veolcity in m/s.

        Returns
        -------
        machNumber : float
            Mach number (unitless).

        """
        
        speedSound = math.sqrt(1.4 * 287 * temp)
        machNumber = velocity / speedSound
        
        return machNumber
    
    ###
    
    def reynoldsNumber(density, velocity, length, mu):
        """
        Calculates Reynolds Number for a component

        Parameters
        ----------
        density : float
            Density at given altitude in kg/m**3.
        velocity : float
            Velocity in m/s.
        length : float
            Average chord length or length of component in m.
        mu : float
            Mu value in (Ns)/m**2.

        Returns
        -------
        reynolds : float
            Reynolds number for the component (unitless).

        """
        reynolds = (density * velocity * length) / mu
        
        return reynolds
    
    ###
    
    def cfcTurbulant(re, mach):
        """
        Calculates the component skin friction coefficient for turbulant airflow.

        Parameters
        ----------
        re : float
            Reynolds Number for the component.
        mach : float
            Mach number.

        Returns
        -------
        cfc : float
            Component skin friction coefficient.

        """
        
        cfc = (0.455) / ( ((math.log10(re)) ** 2.58) * (1 + 0.144 * (mach**2)) ** 0.65 )
        
        return cfc
    
    ###
    
    def ffWing(XCm, tc, mach, angle):
        """
        Calculates the wing form factor.

        Parameters
        ----------
        XCm : float
            Chordwise location of the maximum thickness point of the wing (or similar structure).
            Takes value in m.
        tc : float
            Thickness ratio.
        mach : float
            Mach number.
        angle : float
            Sweepback angle in degrees.

        Returns
        -------
        ffWing : float
            Form factor for the wing (or similar structure) (unitless).

        """
        
        coef1 = (1 + ((0.6 / XCm)*tc) + (100 * ((tc)**4)))
        coef2 = ((1.34 * (mach**0.18)) * (math.cos(math.radians(angle))**0.28))
        ffWing = coef1 * coef2
        
        return ffWing
    
    ###
    
    def sWetAirfoil(tc, sExposed, adjustment = 0):
        """
        Calculates the wetted area of the wing.

        Parameters
        ----------
        tc : float
            Thickness ratio.
        sExposed : float
            Area of the exposed portion of the wing (m**2).
        adjustment : float, optional
            Area of wing as part of the fuselage (m**2).
            Will be removed from sExposed. The default is 0.

        Returns
        -------
        sWetWing : float
            Wetted area of the wing (m**2).

        """
        
        sWetWing = (1.977 + (0.52 * tc))*(sExposed - adjustment)
        
        return sWetWing
    
    ###
    
    def eKnot (AR):
        """
        Calculates the Oswald Efficiency Factor.

        Parameters
        ----------
        AR : float
            Aspect ratio from aircraft wings.

        Returns
        -------
        eKnot : float
            Oswald Efficiency Factor.

        """
        
        eKnot = 1.78 * (1 - (0.045 * (AR**0.68))) - 0.64
        
        return eKnot
    
    ###
    
    def valueK (AR, eKnot):
        """
        Caluclates the induced drag coefficient.

        Parameters
        ----------
        AR : float
            Aspect ratio from aircraft wings.
        eKnot : float
            Oswald Efficiency Factor.

        Returns
        -------
        K : float
            Induced drag coefficient.

        """
        
        K = 1 / (math.pi * AR * eKnot)
        
        return K
    
    ###
    
    def ffFuse(length, aMax):
        """
        Calculates the fuselage (or smooth canopy) form factor.

        Parameters
        ----------
        length : float
            Length of the fuselage (m).
        aMax : float
            Max cross sectional area of the fuselage (m**2).

        Returns
        -------
        ff : float
            Fuselage form factor.

        """
        
        f = length / math.sqrt((math.pi / 4) * aMax)
        ff = 0.9 + (5/f**0.5) + (f/400)
        
        return ff

    ###
    
    def CLCDmax(K, cD0, path):
        """
        Calculates the CL/CD ratio based on specified path.
        If path is set to 1, calculates (CL/CD)max.
        If path is set to 2, calculates (CL(3/2)/CD)max.

        Parameters
        ----------
        K : float
            Induced drag coefficient K.
        cD0 : float
            Total parasite drag.
        path : int
            Function identifier explained above. Changes formula used as necessary

        Returns
        -------
        CLCD : float
            (CL/CD)max value.

        """
        if path == 1:
            CLCD = math.sqrt( 1 / (4 * K * cD0) )
        elif path == 2:
            CLCD = 0.25 * math.pow( (3 / (K * math.pow(cD0, 1/3)) ) , 3/4)
    
        return CLCD
    
    ###
    
    def calcRange(nPR, C, CLCD, W0, W1):
        """
        Calculates the range of an aircraft.

        Parameters
        ----------
        nPR : float
            Prop efficiency ration.
        C : float
            Specific fuel consumption (1 / unit length).
        CLCD : float
            (CL/CD) ratio.
        W0 : float
            Gross takeoff weight.
        W1 : float
            Weight without fuel.

        Returns
        -------
        rangeValue : float
            Range to achieve given CL/CD ratio.

        """
        
        rangeValue = (nPR / C) * CLCD * math.log(W0 / W1)
        
        return rangeValue
    
    ###
    
    def veloToAchieve(rho, W, S, K, cD0, path):
        """
        Calculates the velocity to acheive a cetain optimization based on path number.
        If path is set to 1, calculates velocity for optimal range.
        If path is set to 2, calculates velocity for optimal power.

        Parameters
        ----------
        rho : float
            Air density.
        W : float
            Aircraft weight.
        S : float
            Planform area of the wing.
        K : float
            Induced drag coefficient.
        cD0 : float
            Total parasite drag.
        path : int
            Function identifier explained above. Will change formula if necessary.

        Returns
        -------
        v : float
            Velocity to achieve given parameter.

        """
        
        if path == 1:
            v = math.sqrt( ((2 * W) / (rho * S)) * math.sqrt(K / cD0) )
        elif path == 2:
            v = math.sqrt( ((2 * W) / (rho * S)) * math.sqrt(K / (3 * cD0)) )
        
        return v
    
    ###
    
    def calcEndurance(nPR, rho, S, C, CLCD, W0, W1):
        """
        Calculates the endurance for a given CL/CD ratio

        Parameters
        ----------
        nPR : float
            Prop efficiency ratio.
        rho : float
            Air density.
        S : float
            Wing planform area.
        C : float
            Specific fuel consumption (1 / unit length).
        CLCD : float
            CL/CD ratio.
        W0 : float
            Gross takeoff weight.
        W1 : TYPE
            No fuel weight.

        Returns
        -------
        Endurance : float
            Endurance of an aircraft for given CL/CD ratio.

        """
        
        Endurance = ((nPR * math.sqrt(2 * rho * S)) / C) * CLCD * ( (1/math.sqrt(W1)) - (1/math.sqrt(W0)) )
    
        return Endurance
    

