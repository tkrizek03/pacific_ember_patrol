import math
import numpy as np
import matplotlib.pyplot as plt

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
    
    ###
    
    def liftCalc(rho, v, sWing, CL):
        
        lift = .5 * rho * (v**2) * sWing * CL
        
        return lift
    
    ###
    
    def dragCalc(rho, v, sWing, cD0, K, CL):
        
        drag = 0.5 * rho * (v**2) * sWing * (cD0 + (K*(CL**2)))
        
        return drag
    
    ###
    
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
    
W0 = 85143.5 # N
sWing = 44.3 # m**2
CLmax = 1.383 
K = 0.0594
cD0 = 0.018
nPR = 0.9
totalP = 932000*2 # W
LDmax = 15.2
vCruise = 89.408 # m/s

pASL = nPR * totalP
print("Power Available at Sea Level: "f"{pASL}"" [W]")

# Initializes empty lists
tempList = []
densityList = []
pressureList = []
vStallList = []
vRCMaxList = []
pRCruiseList = []
pRList = []
pAList = []
pEList =[]
pEOptList = []
RCList = []
RCMaxList = []

y = np.arange(25001) # Creates array from 1 to 25000


for i in range(25001):
    # Loops from 0 m to 25000 m
    
    # Determines temp for each loop iteration, appends to list
    t = AircraftFormulas.StdAtm(i)[1]
    tempList.append(t)
    
    # Determines density for each loop iteration, appends to list
    d = AircraftFormulas.StdAtm(i)[2]
    densityList.append(d)
    
    # Determines pressure for each loop iteration, appends to list
    p = AircraftFormulas.StdAtm(i)[3]
    pressureList.append(p)
    
    # Determines stall velocity for each loop iteration, appends to list
    vStall = math.sqrt((2 * W0) / (d * sWing * CLmax))
    vStallList.append(vStall)
    
    # Determines max R/C velocity for each loop iteration, appends to list
    vRCMax = math.sqrt( (2/d) * math.sqrt( K / (3*cD0) ) * (W0/sWing) )
    vRCMaxList.append(vRCMax)
    
    # Determines power required for cruise velocity for each loop iteration, appends to list
    pRCruise = ( .5 * d * (vCruise**3) * sWing * cD0 ) + ( (2 * K * (W0**2)) / (d * vCruise * sWing) )
    pRCruiseList.append(pRCruise)
    
    # Determines min power required for each loop iteration, appends to list
    pR = ( .5 * d * (vRCMax**3) * sWing * cD0 ) + ( (2 * K * (W0**2)) / (d * vRCMax * sWing) )
    pRList.append(pR)
    
    # Determines power available for each loop iteration, appends to list
    pA = pASL * (d/1.225)
    pAList.append(pA)
    
    # Determines excess power for each loop iteration, appends to list
    pE = pA - pRCruise
    pEList.append(pE)
    
    # Determines max excess power for each loop iteration, appends to list
    pEOptimal = pA - pR
    pEOptList.append(pEOptimal)
    
    # Determines R/C for each loop iteration, appends to list
    RC = pE / W0
    RCList.append(RC)
    
    # Determines R/C max for each loop iteration, appends to list
    RCMax = pEOptimal / W0
    RCMaxList.append(RCMax)


plt.figure(figsize=(6,6))
plt.title('Power Relationships at Altitude')
plt.ylabel('Power [W]')
plt.xlabel('Altitude [m]')
plt.plot(y, pRCruiseList, 'b', label="Power Required at Cruise")
plt.plot(y, pAList, 'g', label="Power Available")
plt.plot(y, pEList, 'r', label="Excess Power")
plt.xlim([0, 15000])
plt.ylim([0, 1500000])
plt.ticklabel_format(useOffset=False, style='plain')
plt.legend()
plt.show()

# Creates new R/C max list with only values above 0
newRCMaxList = [value for value in RCMaxList if value > 0]
newRCMaxList.append(0) # Adds final zero value to end of list
numPoints = len(newRCMaxList) # Finds number of R/C values above zero 
newY = np.arange(numPoints) # Creates new array for the length of the new R/C list

targetValue = (100 / (3.2808399*60)) # Creates target service ceiling value in terms of m/s
serviceCeiling = min(newRCMaxList, key=lambda x: abs(x - targetValue)) # Finds the value in new RC list closest to service ceiling
serviceCeilingIndex = newRCMaxList.index(serviceCeiling) # Finds location of service ceiling in new RC list
serviceCeilingAltitude = serviceCeilingIndex + 1 # Accounts for index values starting at 0
print("Service Ceiling Altitude: "f"{serviceCeilingAltitude}"" [m]")

absoluteCeilingAltitude = len(newRCMaxList) # Finds length of new RC list as absolute ceiling
print("Absolute Ceiling Altitude: "f"{absoluteCeilingAltitude}"" [m]")

# RC plot with ceiling lines
plt.figure(figsize=(6,6))
plt.title('R/C vs. Altitude')
plt.xlabel('Altitude [m]')
plt.ylabel('R/C [m/s]')
plt.plot(newY, newRCMaxList, label="R/C Max")
plt.axvline(x=serviceCeilingIndex + 1, color='g', label="Service Ceiling")
plt.axvline(x=len(newRCMaxList), color='r', label="Absolute Ceiling")
plt.xlim([0, 10000])
plt.ylim([0, 20])
plt.legend()
plt.show()

#############################

"""
Pacific Ember Patrol

Calcualtion interface for AircraftFormulas class.
"""
################
# Variable definitions and calculations for atmospheric conditions and base wing/weight values.
t = 272.3167 # Temp in kelvin
rho = 0.9629 # Density in kg/m**3
v = 89.408 # Cruising velocity in m/s

Wp = 1000 * 8.8707 * 9.81 / 2.205  # Payload weight in N (capacity * component weight * conversion)
We = 7000 * 9.81 / 2.205 # Empty weight in N 
Wf = 6.8 * 480.45 * 9.81  / 2.205 # Fuel weight in N (capacity * component weight * conversion)
W0 = Wp + We + Wf # Total weight in N
fracWp = Wp/W0
fracWe = We/W0
fracWf = Wf/W0
sWing = 44.3275
print("Empty Weight: " f"{We:.6}"" [N]")
print("Fuel Weight: " f"{Wf:.6}"" [N]")
print("Payload Weight: " f"{Wp:.6}"" [N]")
print("Total Weight: " f"{W0:.6}"" [N]")
print("Empty Weight Fraction: " f"{fracWe:.2}")
print("Payload Weight Fraction: " f"{fracWp:.2}")
print("Fuel Weight Fraction: " f"{fracWf:.2}")

cAvg = 2.67 # Average chord length in m
b = sWing / cAvg # Wingspan in m
AR = (b**2) / sWing # Aspect ratio (unitless)
print("Wing Planform Area: " f"{sWing:.3}"" [m**2]")
print("Average Chord Length: " f"{cAvg:.3}"" [m]")
print("Wingspan: " f"{b:.4}"" [m]")
print("Aspect Ratio: " f"{AR:.3}")

XCm = (cAvg * 0.309) # Chordwise location of maximum thickness point in m
tc = 0.15 # Thickness ratio
angle = 0 # Sweepback angle
sExposed = sWing # Exposed area of the wing (m**2)
qWing = 1.0 # Interference factor of the wing
print("Chordwise Location of Max Thickness: "f"{XCm:.3}"" m")
print("Sweepback Angle: "f"{angle}"" degrees")


################
# Performs calculates for variables based on cruising conditions

mu = AircraftFormulas.dynamicViscosity(t) # Mu value
print("Mu Value: " f"{mu:.2e}"" [(Ns)/m**2]")

reWing = AircraftFormulas.reynoldsNumber(rho, v, cAvg, mu) # Wing Reynolds Number
print("Wing Reynolds Number: "f"{reWing:.2e}")

m = AircraftFormulas.mach(t, v) # Mach number
print("Mach Number: "f"{m:.3}")

eKnot =  AircraftFormulas.eKnot(AR)
print("Oswald Efficiency Factor: "f"{eKnot:.3}")
K = AircraftFormulas.valueK(AR, eKnot)
print("Drag Due to Lift Value (K): "f"{K:.3}")


################
# Performs calculations for the wing drag coefficient

cfcWing = AircraftFormulas.cfcTurbulant(reWing, m)
print("Wing Skin Friction Coefficient: "f"{cfcWing:.2}")
ffWing = AircraftFormulas.ffWing(XCm, tc, m, angle)
print("Wing Form Factor: "f"{ffWing:.3}")
sWetWing = AircraftFormulas.sWetAirfoil(tc, sExposed)
print("Wetted Wing Area: "f"{sWetWing:.3}"" [m**2]")

cD0wing = AircraftFormulas.coefficientDragComponent(cfcWing, ffWing, qWing, sWetWing, sWing)
print("Wing Drag Coefficient: "f"{cD0wing:.2}")

################
# 2D/3D lift curve slope adjustments and values (updated as needed)

a = 11.83 # Placeholder. To be updated as needed. Entered as degrees.
aL0 = -4.17
a2D = 0.1224
a3D = a2D / (1 + ( (57.3*a2D) / (math.pi*eKnot*AR) ) )
Cl = a2D*(a - aL0)
CL = a3D*(a - aL0)
print("2D Lift Curve Max CL: "f"{Cl:.4}")
print("3D Lift Curve Max CL: "f"{CL:.4}")

################
# CL cruise

CLcruise = W0 / (.5 * rho * (v**2) * sWing)
print("CL for Cruise Conditions: "f"{CLcruise:.4}")

################
# CD,0,fuselage calculations

length = 11
aMax = 3.75
aTop = 10.77
aSide = 9.19
sWetFuse = 3.4 * ((aTop + aSide)/2)
qFuse = 1

reFuse = AircraftFormulas.reynoldsNumber(rho, v, length, mu)

cfcFuse = AircraftFormulas.cfcTurbulant(reFuse, m)
ffFuse = AircraftFormulas.ffFuse(length, aMax)

cD0fuse = AircraftFormulas.coefficientDragComponent(cfcFuse, ffFuse, qFuse, sWetFuse, sWing)

print("Fuselage Reynolds Number: "f"{reFuse:.2e}")
print("Fuselage Skin Friction Coefficient: "f"{cfcFuse:.2}")
print("Fuselage Form Factor: "f"{ffFuse:.3}")
print("Wetted Fuselage Area: "f"{sWetFuse:.3}"" [m**2]")

print("Fuselage Drag Coefficient: "f"{cD0fuse:.2}")

################
# CD,0,stabilizer calculations

cAvgStabilizer = 1.34
XCmStabilizer = 0.268
tcStabilizer = .2
sHStabilizer = 4.25 
sVStabilizer = 3.1
sHStabilizerTotal = sHStabilizer * 2
angleStabilizer = 0 
qStabilizer = 1.05

reStabilizer = AircraftFormulas.reynoldsNumber(rho, v, cAvgStabilizer, mu)

sWetHStabilizer = AircraftFormulas.sWetAirfoil(tcStabilizer, sHStabilizerTotal)
sWetVStabilizer = AircraftFormulas.sWetAirfoil(tcStabilizer, sVStabilizer)
cfcStabilizer = AircraftFormulas.cfcTurbulant(reStabilizer, m)
ffStabilizer = AircraftFormulas.ffWing(XCmStabilizer, tcStabilizer, m, angleStabilizer)

cD0HStabilizer = AircraftFormulas.coefficientDragComponent(cfcStabilizer, ffStabilizer, qStabilizer, sWetHStabilizer, sWing)
cD0VStabilizer = AircraftFormulas.coefficientDragComponent(cfcStabilizer, ffStabilizer, qStabilizer, sWetVStabilizer, sWing)


print("Empennage Reynolds Number: "f"{reStabilizer:.2e}")
print("Empennage Skin Friction Coefficient: "f"{cfcStabilizer:.2}")
print("Empennage Form Factor: "f"{ffStabilizer:.3}")
print("Horizontal Stabilizer Wetted Area: "f"{sWetHStabilizer:.3}"" [m**2]")
print("Vertical Stabilizer Wetted Area: "f"{sWetVStabilizer:.3}"" [m**2]")

print("Horizontal Stabilizer Drag Coefficient: "f"{cD0HStabilizer:.2}")
print("Vertical Stabilizer Drag Coefficient: "f"{cD0VStabilizer:.2}")

################
# CD,0,landing gear calculations

cDGear = 0.25
sFrontalMain = 229.5 / 1550
sFrontalRear = 108 / 1550
cDStrut = 0.3
sFrontalStrut = 10.5 / 1550
qGear = 1.3
qStrut = 1

cD0MainGear = cDGear * qGear * (sFrontalMain/sWing)
cD0RearGear = cDGear * qGear * (sFrontalRear/sWing)
cD0Struts = cDStrut * qStrut * (sFrontalStrut/sWing)
totalcD0gear = (2 * cD0MainGear) + (cD0RearGear) + (3 * cD0Struts)

print("Main Gear Drag Coefficient: "f"{cD0MainGear:.2}")
print("Rear Gear Drag Coefficient: "f"{cD0RearGear:.2}")
print("Strut Drag Coefficient: "f"{cD0Struts:.2}")

print("Total Gear Drag Coefficient: "f"{totalcD0gear:.2}")

################
# CD,0,nacelle calculations

aMaxN = .28
sWetN = 2.54
qN = 1.5
lenN = 1.9

reN = AircraftFormulas.reynoldsNumber(rho, v, lenN, mu)

cfcN = AircraftFormulas.cfcTurbulant(reN, m)
fN = lenN / (math.sqrt((4 * aMaxN) / math.pi))
ffN = 1 + (0.35 / fN)

cD0Nacelle = 2 * AircraftFormulas.coefficientDragComponent(cfcN, ffN, qN, sWetN, sWing)

print("Total Nacelle Drag Coefficient: "f"{cD0Nacelle:.2}")

################
# CD,0,total calculations

cD0total = cD0wing + cD0fuse + cD0HStabilizer + cD0VStabilizer + totalcD0gear + cD0Nacelle

print("Total Parasite Drag Coefficient: "f"{cD0total:.2}")

################
# Range calculations (assumes no payload is dropped)
vStall = math.sqrt((2 * W0) / (rho * sWing * CL))

W1 = We + Wp
SFC = 0.54 # lb/hp*hr
C = (SFC * 608.22739 * 9.81*2) / (1000 * 1000 * 3600) 
nPR = 0.9
CLCDmax = AircraftFormulas.CLCDmax(K, cD0total, 1)

vRFull = AircraftFormulas.veloToAchieve(rho, W0, sWing, K, cD0total, 1)
vREmpty = AircraftFormulas.veloToAchieve(rho, W1, sWing, K, cD0total, 1)

maxRangeM = AircraftFormulas.calcRange(nPR, C, CLCDmax, W0, W1)
maxRangeKM = maxRangeM / 1000

print("Brake Specific Fuel Consumption: "f"{C:.2e}"" [1/m]")
print("Empty Fuel Weight: "f"{W1:.6}"" [N]")
print("(CL/CD)max: "f"{CLCDmax:.3}")
print("Stall Velocity: "f"{vStall:.4}"" [m/s]")
print("Velocity to Achieve Max Range with Fuel: "f"{vRFull:.4}"" [m/2]")
print("Velocity to Achieve Max Range with No Fuel: "f"{vREmpty:.4}"" [m/s]")
print("Maximum Range: "f"{maxRangeKM:.5}"" [km]")


################
# Power calculations (assumes no payload is dropped)
CL32CDmax = AircraftFormulas.CLCDmax(K, cD0total, 2)

vEFull = AircraftFormulas.veloToAchieve(rho, W0, sWing, K, cD0total, 2)
vEEmpty = AircraftFormulas.veloToAchieve(rho, W1, sWing, K, cD0total, 2)

maxEnduranceS = AircraftFormulas.calcEndurance(nPR, rho, sWing, C, CL32CDmax, W0, W1)
maxEnduranceHR = maxEnduranceS / 3600

print("(CL(3/2)/CD)max: "f"{CL32CDmax:.3}")
print("Velocity to Achieve Max Endurance with Fuel: "f"{vEFull:.4}"" [m/2]")
print("Velocity to Achieve Max Endurance with No Fuel: "f"{vEEmpty:.4}"" [m/s]")
print("Maximum Endurance: "f"{maxEnduranceHR:.3}"" [hr]")

################
# Begin Part F Calculations #
# Pull-Up Maneuver Calculations
g = 9.81 
nPosStruc = 3.8
nPosAero = (.5 * rho * (v**2) * sWing * CL) / W0
print("Structural Limit: "f"{nPosStruc:.2}")
print("Aerodynamic Limit: "f"{nPosAero:.3}")

rPullUpS = (v**2) / (g * (nPosStruc - 1))
rPullUpA = (v**2) / (g * (nPosAero - 1))
print("Structurally Limited Pull Up Radius: "f"{rPullUpS:.5}"" [m]")
print("Aerodynamically Limited Pull Up Radius: "f"{rPullUpA:.5}"" [m]")

pullUpRate = (g * (nPosAero - 1)) / v
print("Aerodynamically Limited Pull Up Rate: "f"{pullUpRate:.2}"" [radians/s]")

################
# Level Turn Maneuver Calculations
rTurnS = (v**2) / (g * math.sqrt((nPosStruc)**2 - 1))
rTurnA = (v**2) / (g * math.sqrt((nPosAero)**2 - 1))
print("Structurally Limited Level Turn Radius: "f"{rTurnS:.5}"" [m]")
print("Aerodynamically Limited Level Turn Radius: "f"{rTurnA:.5}"" [m]")

levelTurnRate = (g * math.sqrt((nPosAero**2)-1) / (v))
print("Aerodynamically Limited Level Turn Rate: "f"{levelTurnRate:.2}"" [radians/s]")

################
# Maneuvering Velocity Calculation
vA = math.sqrt( (2*nPosStruc*W0) / (rho*CL*sWing))
print("Maximum Velocity for Maneuvering: "f"{vA:.5}"" [m/s]")

################
# Takeoff Calculation
totalP = 932000*2
pASL = totalP * nPR
muR = 0.03
rhoT = 1.0984
vLO = 1.2*vStall
aRolling = 3
aClimb = 3
nTakeoff = 1.15
CLrolling = a3D*(aRolling - aL0)
liftT = AircraftFormulas.liftCalc(rhoT, 0.7*vLO, sWing, CLrolling)
dragT = AircraftFormulas.dragCalc(rhoT, 0.7*vLO, sWing, cD0total, K, CLrolling)
# pRTakeoff = ( .5 * rhoT * ((0.7*vLO)**3) * sWing * cD0total ) + ( (2 * K * (W0**2)) / (rhoT * (0.7*vLO) * sWing) )
pAT = pASL * (rhoT/1.225)
thrustT = pAT / (0.7*vLO)
print("Coefficient of Lift while Rolling: "f"{CLrolling:.3}")

takeoffGroundRoll = (1.44 * (W0**2)) / (g * rhoT * sWing * CL * (thrustT - dragT - (muR * (W0 - liftT))))
print("Lift Generated During Takeoff Ground Roll: "f"{liftT:.6}"" [N]")
print("Drag Generated During Takeoff Ground Roll: "f"{dragT:.5}"" [N]")
print("Takeoff Ground Roll Distance: "f"{takeoffGroundRoll:.5}"" [m]")

rTakeoff = (vLO**2) / (g * (nTakeoff - 1))
takeoffTransition = rTakeoff * math.sin((math.pi*aClimb)/180)
print("Takeoff Maneuver Radius: "f"{rTakeoff:.5}"" [m]")
print("Takeoff Transition Distance: "f"{takeoffTransition:.5}"" [m]")

hATakeoff = (35/3.281) - rTakeoff + (rTakeoff * math.cos((math.pi*aClimb)/180))
takeoffAir = hATakeoff / (math.tan((math.pi*aClimb)/180))

print("Takeoff Air Distance: "f"{takeoffAir:.5}"" [m]")

totalTakeoff = takeoffGroundRoll + takeoffTransition + takeoffAir
print("Total Takeoff Distance: "f"{totalTakeoff:.5}"" [m]")

################
# Landing Calculation
vTD = 1.3*vStall
vFlare = vTD
thetaFlare = 3
thetaApproach = thetaFlare
muB = 0.5
nLanding = 1.12 # From difference of 0.03 in HW9

rLanding = (vFlare**2) / (g * (nLanding - 1))
hALanding = (50/3.281) - rLanding + (rLanding * (math.cos(math.pi*thetaApproach/180)))
landingAir = hALanding / math.tan(math.pi*thetaFlare/180)
print("Landing Maneuver Radius: "f"{rLanding:.5}"" [m]")
print("Landing Air Distance: "f"{landingAir:.5}"" [m]")

landingFlare = rLanding * math.sin(math.pi*thetaFlare/180)
print("Landing Flare Distance: "f"{landingFlare:.5}"" [m]")

liftL = AircraftFormulas.liftCalc(rhoT, 0.7*vTD, sWing, CLrolling)
dragL = AircraftFormulas.dragCalc(rhoT, 0.7*vTD, sWing, cD0total, K, CLrolling)
landingGroundRoll = (1.69 * (W1**2)) / (g * rhoT * sWing * CL * (dragL + (muB * (W1 - liftL))))
print("Lift Generated During Landing Ground Roll: "f"{liftL:.6}"" [N]")
print("Drag Generated During Landing Ground Roll: "f"{dragL:.5}"" [N]")
print("Landing Ground Roll Distance: "f"{landingGroundRoll:.5}"" [m]")

totalLanding = landingAir + landingFlare + landingGroundRoll
print("Total Landing Distance: "f"{totalLanding:.5}"" [m]")

################
# Fuel Consumption Calculations
# Landing

W1Takeoff = W0 / (math.exp( ( totalTakeoff * C * dragT ) / (nPR * liftT)))
takeoffFuelConsumption = ((W0 - W1Takeoff) * 2.205) / (6.8 * 9.81)
print("Fuel Consumed During Takeoff: "f"{takeoffFuelConsumption:.3}"" [gal]")

# Climb (AltDependantCalc)
climbRange = 30 * 1.609 * 1000
minClimbAlt = 1122
maxClimbAlt = 10000 * .3048
avgClimbAlt = int(round((minClimbAlt + maxClimbAlt) / 2))
liftC = W0 * math.cos(math.radians(3))
avgPR = pRList[avgClimbAlt]
avgDensity = densityList[avgClimbAlt]
pAC = pASL * (avgDensity/1.225)
powerCFraction = avgPR / pAC
CClimb = powerCFraction * C
vClimb = RCList[avgClimbAlt] / math.sin(math.radians(3))
thrustC = pAC / vClimb
dragC = thrustC - (W0 * math.sin(math.radians(3)))

W1Climb = W0 / (math.exp( ( climbRange * CClimb * dragC ) / (nPR * liftC)))
climbFuelConsumption = ((W0 - W1Climb) * 2.205) / (6.8 * 9.81)
print("Range of Climb: "f"{climbRange:.6}"" [m]")
print("Fuel Consumed During Climb: "f"{climbFuelConsumption:.3}"" [gal]")

# Cruise
# Assumes loiter flight of 90 miles
cruiseRange = 90 * 1.609 * 1000
powerCrFraction = pRList[3048] / pAList[3048]
CCruise = powerCrFraction * C
CDcruise = cD0total + (K * (CLcruise**2))
dragCruise = CDcruise * .5 * densityList[3048] * (v**2) * sWing
liftCruise = CLcruise * .5 * densityList[3048] * (v**2) * sWing

W1Cruise = W0 / (math.exp( ( cruiseRange * CCruise * dragCruise) / (nPR * liftCruise)))
cruiseFuelConsumption = ((W0 - W1Cruise) * 2.205) / (6.8 * 9.81)
print("Fuel Consumed During Loiter @ Cruise Conditions: "f"{cruiseFuelConsumption:.3}"" [gal]")

