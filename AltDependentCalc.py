import math
import numpy as np
import matplotlib.pyplot as plt

import StandAtmFunc

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
    t = StandAtmFunc.StdAtm(i)[1]
    tempList.append(t)
    
    # Determines density for each loop iteration, appends to list
    d = StandAtmFunc.StdAtm(i)[2]
    densityList.append(d)
    
    # Determines pressure for each loop iteration, appends to list
    p = StandAtmFunc.StdAtm(i)[3]
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

plt.figure(figsize=(6,6))
plt.title('R/C vs. Altitude')
plt.xlabel('Altitude [m]')
plt.ylabel('R/C [m/s]')
plt.plot(y, RCList)
plt.xlim([0, 10000])
plt.ylim([0, 20])
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

rhoFinder = densityList[1123]
