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
Wf = 6.8 * 480.45 * 9.81  / 2.205 # Fuel weight in kg (capacity * component weight * conversion)
W0 = Wp + We + Wf # Total weight in kg
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
# CD,0,total calculations

cD0total = cD0wing + cD0fuse + cD0HStabilizer + cD0VStabilizer + totalcD0gear

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
print("Maneuvering Velocity: "f"{vA:.5}"" [m/s]")

################
# Takeoff Calculation
totalP = 932000*2
pASL = totalP * nPR
muR = 0.03
rhoT = 1.0984
vLO = 1.2*vStall
aRolling = 6
aClimb = 3
nTakeoff = 1.15
CLrolling = W0 / (.5 * rhoT * (vLO**2) * sWing)
# CLrolling = a3D*(aRolling - aL0)
# print(CLrolling)
liftT = AircraftFormulas.liftCalc(rhoT, 0.7*vLO, sWing, CLrolling)
dragT = AircraftFormulas.dragCalc(rhoT, 0.7*vLO, sWing, cD0total, K, CLrolling)
pAT = pASL * (rhoT/1.225)
thrustT = pAT / (0.7*vLO)

takeoffGroundRoll = (1.44 * (W0**2)) / (g * rhoT * sWing * CL * (thrustT - dragT - (muR * (W0 - liftT))))
print("Takeoff Ground Roll Distance: "f"{takeoffGroundRoll:.5}"" [m]")

rTakeoff = (vLO**2) / (g * (nTakeoff - 1))
takeoffTransition = rTakeoff * math.sin((math.pi*aClimb)/180)
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
print("Landing Air Distance: "f"{landingAir:.5}"" [m]")

landingFlare = rLanding * math.sin(math.pi*thetaFlare/180)
print("Landing Flare Distance: "f"{landingFlare:.5}"" [m]")

liftL = AircraftFormulas.liftCalc(rhoT, 0.7*vTD, sWing, CLrolling)
dragL = AircraftFormulas.dragCalc(rhoT, 0.7*vTD, sWing, cD0total, K, CLrolling)
landingGroundRoll = (1.69 * (We**2)) / (g * rhoT * sWing * CL * (dragL + (muB * (We - liftL))))
print("Landing Ground Roll Distance: "f"{landingGroundRoll:.5}"" [m]")

totalLanding = landingAir + landingFlare + landingGroundRoll
print("Total Landing Distance: "f"{totalLanding:.5}"" [m]")




