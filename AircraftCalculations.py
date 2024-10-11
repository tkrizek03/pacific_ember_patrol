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

Wp = 1000 * 8.8707 * 9.81 / 2.205  # Payload weight in kg (capacity * component weight * conversion)
We = 7000 * 9.81 / 2.205 # Empty weight in kg 
Wf = 6.8 * 450 * 9.81  / 2.205 # Fuel weight in kg (capacity * component weight * conversion)
W0 = Wp + We + Wf # Total weight in kg
fracWp = Wp/W0
fracWe = We/W0
fracWf = Wf/W0
sWing = W0 / 1900
print("Empty Weight: " f"{We:.6}"" N")
print("Fuel Weight: " f"{Wf:.6}"" N")
print("Payload Weight: " f"{Wp:.6}"" N")
print("Total Weight: " f"{W0:.6}"" N")
print("Empty Weight Fraction: " f"{fracWe:.2}")
print("Payload Weight Fraction: " f"{fracWp:.2}")
print("Fuel Weight Fraction: " f"{fracWf:.2}")

cAvg = 2.67 # Average chord length in m
b = sWing / cAvg # Wingspan in m
AR = (b**2) / sWing # Aspect ratio (unitless)
print("Wing Planform Area: " f"{sWing:.3}"" m**2")
print("Average Chord Length: " f"{cAvg:.3}"" m")
print("Wingspan: " f"{b:.4}"" m")
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
print("Drag Due to Lift Value (K): "f"{K:.3}")


################
# Performs calculations for the wing drag coefficient

cfcWing = AircraftFormulas.cfcTurbulant(reWing, m)
print("Wing Skin Friction Coefficient: "f"{cfcWing:.2}")
ffWing = AircraftFormulas.ffWing(XCm, tc, m, angle)
print("Wing Form Factor: "f"{ffWing:.3}")
sWetWing = AircraftFormulas.sWetAirfoil(tc, sExposed)
print("Wetted Wing Area: "f"{sWetWing:.3}"" m**2")

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
aMax = 4.32
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
print("Wetted Fuselage Area: "f"{sWetFuse:.3}"" m**2")

print("Fuselage Drag Coefficient: "f"{cD0fuse:.2}")

################
# CD,0,stabilizer calculations

cAvgStabilizer = 1.34
XCmStabilizer = 0.268
tcStabilizer = .2
sStabilizer = 4.25 
sStabilizerTotal = sStabilizer * 3
angleStabilizer = 0 
qStabilizer = 1.05

reStabilizer = AircraftFormulas.reynoldsNumber(rho, v, cAvgStabilizer, mu)

sWetStabilizer = AircraftFormulas.sWetAirfoil(tcStabilizer, sStabilizerTotal)
cfcStabilizer = AircraftFormulas.cfcTurbulant(reStabilizer, m)
ffStabilizer = AircraftFormulas.ffWing(XCmStabilizer, tcStabilizer, m, angleStabilizer)


cD0Stabilizer = AircraftFormulas.coefficientDragComponent(cfcStabilizer, ffStabilizer, qStabilizer, sStabilizerTotal, sWing)

print("Empennage Reynolds Number: "f"{reStabilizer:.2e}")
print("Empennage Skin Friction Coefficient: "f"{cfcStabilizer:.2}")
print("Empennage Form Factor: "f"{ffStabilizer:.3}")
print("Total Wetted Empennage Area: "f"{sWetStabilizer:.3}"" m**2")

print("Empennage Drag Coefficient: "f"{cD0Stabilizer:.2}")


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
# CD,0,total calculations

cD0total = cD0wing + cD0fuse + cD0Stabilizer + totalcD0gear

print("Total Parasite Drag Coefficient: "f"{cD0total:.2}")


