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
