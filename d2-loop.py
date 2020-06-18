#!/usr/bin/python3

import CoolProp.CoolProp as CP
import numpy as np
fluid='Deuterium'

T=20 # K
P=101325 # Pa
h=CP.PropsSI('H','T',T,'P',P,fluid)
s=CP.PropsSI('S','T',T,'P',P,fluid)
print(T,P,h,s)

T=300 # K
P=101325 # Pa
h=CP.PropsSI('H','T',T,'P',P,fluid)
s=CP.PropsSI('S','T',T,'P',P,fluid)
print(T,P,h,s)

P=101325 # Pa
T=CP.PropsSI('T','P',P,'Q',1,fluid)
h=CP.PropsSI('H','P',P,'Q',1,fluid)
s=CP.PropsSI('S','P',P,'Q',1,fluid)
# Q=1=vapour
print(T,P,h,s)

P=101325 # Pa
T=CP.PropsSI('T','P',P,'Q',0,fluid)
h=CP.PropsSI('H','P',P,'Q',0,fluid)
s=CP.PropsSI('S','P',P,'Q',0,fluid)
# Q=0=liquid
print(T,P,h,s)
# conclusion: they set the enthalpy and entropy to zero for liquid at
# saturation at standard pressure, by default, for this fluid.  This
# seems to be consistent with one of the conventions in CoolProp.

for T in np.arange(18.624,25,0.01):
    d=CP.PropsSI('D','T',T,'Q',0,fluid)
    p=CP.PropsSI('P','T',T,'Q',0,fluid)
    print(T,d,p)

# Based on Cam's design
    
p_psi=20. # PSI
p=p_psi*6894.76 # Pa

# What's the temperature at which the two-phase gas-liquid condition exists (saturated state)

Tmax=CP.PropsSI('T','P',p,'Q',0,fluid)
print(Tmax)

# What's the temperature at which it solidifies?

#p=101325

state=CP.AbstractState("","Deuterium")
print(state.has_melting_line())
# Conclusion:  seems like CoolProp does not have this information.

p_psi=20. # PSI
p=p_psi*6894.76 # Pa
T=23.0 # K

h1=CP.PropsSI('H','T',T,'P',p,fluid)
phase=CP.PhaseSI('T',T,'P',p,fluid)
print(T,p,h1,phase)

T=24.0 # K
h2=CP.PropsSI('H','T',T,'P',p,fluid)
density=CP.PropsSI('D','T',T,'P',p,fluid)
phase=CP.PhaseSI('T',T,'P',p,fluid)
print(T,p,h2,phase)

qdot=50 # W
delta_t=0.125*density*(h2-h1)/qdot
print(delta_t/60)

mdot=qdot/(h2-h1) # kg/s
print(mdot)


T=20.5 # K
p_psi=20. # PSI
p=p_psi*6894.76 # Pa
viscosity=CP.PropsSI('VISCOSITY','T',T,'P',p,'hydrogen')
print(viscosity)
