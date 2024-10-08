import math
from AircraftFormulas import *

"""
Pacific Ember Patrol

Calcualtion interface for AircraftFormulas class.
"""
################
# Variable definitions and calculations for atmospheric conditions and base wing/weight values.
t = 272.3167 # Temp in kelvin
rho = 0.9629 # Density in kg/m**3
v = 89.408 # Cruising velocity in m/s

Wp = 1000 * 8.8707 / 2.205 # Payload weight in kg (capacity * component weight * conversion)
We = 7000 / 2.205 # Empty weight in kg 
Wf = 6.8 * 450 / 2.205 # Fuel weight in kg (capacity * component weight * conversion)
W0 = Wp + We + Wf # Total weight in kg
fracWp = Wp/W0
fracWe = We/W0
fracWf = Wf/W0
weightForce = W0 * 9.81
sWing = weightForce / 1900
print("Empty Weight: " f"{We:.6}"" kg")
print("Fuel Weight: " f"{Wf:.6}"" kg")
print("Payload Weight: " f"{Wp:.6}"" kg")
print("Total Weight: " f"{W0:.6}"" kg")
print("Empty Weight Fraction: " f"{fracWe:.2}")
print("Payload Weight Fraction: " f"{fracWp:.2}")
print("Fuel Weight Fraction: " f"{fracWf:.2}")

cAvg = 2.67 # Average chord length in m
b = sWing / cAvg # Wingspan in m
AR = (b**2) / sWing # Aspect ratio (unitless)
print("Wing Planform Area: " f"{sWing:.3}"" m**2")
print("Average Chord Length: " f"{cAvg:.3}"" m")
print("Wingspan: " f"{b:.2e}"" m")
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
print("Mu Value: " f"{mu:.2e}"" (Ns)/m**2")

reWing = AircraftFormulas.reynoldsNumber(rho, v, cAvg, mu) # Wing Reynolds Number
print("Wing Reynolds Number: "f"{reWing:.2e}")

m = AircraftFormulas.mach(t, v) # Mach number
print("Mach Number: "f"{m:.3}")

eKnot =  AircraftFormulas.eKnot(AR)
print("Oswald Efficiency Factor: "f"{eKnot:.3}")
K = AircraftFormulas.valueK(AR, eKnot)
print("Drag Due to Lift Value (K): "f"{K:.2e}")


################
# Performs calculations for the wing drag coefficient

cfcWing = AircraftFormulas.cfcTurbulant(reWing, m)
print("Wing Skin Friction Coefficient: "f"{cfcWing:.2e}")
ffWing = AircraftFormulas.ffWing(XCm, tc, m, angle)
print("Wing Form Factor: "f"{ffWing:.3}")
sWetWing = AircraftFormulas.sWetWing(tc, sExposed)
print("Wetted Wing Area: "f"{sWetWing:.4}"" m**2")

cD0wing = AircraftFormulas.coefficientDragComponent(cfcWing, ffWing, qWing, sWetWing, sWing)
print("Wing Drag Coefficient: "f"{cD0wing:.2e}")

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

CLcruise = weightForce / (.5 * rho * (v**2) * sWing)
print("CL for Cruise Conditions: "f"{CLcruise:.4}")

################
# CD,0,fuselage calculations

length = 0 
aMax = 0
sWetFuse = 0
qFuse = 1

reFuse = AircraftFormulas.reynoldsNumber(rho, v, length, mu)

cfcFuse = AircraftFormulas.cfcTurbulant(reFuse, m)
ffFuse = AircraftFormulas.ffFuse(length, aMax)

cD0fuse = AircraftFormulas.coefficientDragComponent(cfcFuse, ffFuse, qFuse, sWetFuse, sWing)

################
# CD,0,horizontal stabilizer calculations

cAvgHorizontal = 0
XCmHorizontal = 0
tcHorizontal = 0
sHorizontal = 0 
sExposedHorizontal = sHorizontal
angleHorizontal = 0 
qHorizontal = 0

reHorizontal = AircraftFormulas.reynoldsNumber(rho, v, cAvgHorizontal, mu)

cfcHorizontal = AircraftFormulas.cfcTurbulant(reHorizontal, m)
ffHorizontal = AircraftFormulas.ffWing(XCmHorizontal, tcHorizontal, m, angleHorizontal)
sWetHorizontal = AircraftFormulas.sWetWing(tcHorizontal, sExposedHorizontal)

cD0horizontal = AircraftFormulas.coefficientDragComponent(cfcHorizontal, ffHorizontal, qHoriztontal, sWetHorizontal, sWing)