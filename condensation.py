#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import numpy as np
from scipy import interpolate
from scipy.integrate import odeint

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
fluid='Deuterium'
g=9.8 #m/s

beta_t=.012 # K^(-1), relative slope of density with temperature (a la
# Boussinesq)


d2=4.76*0.0254 # (m) inner diameter of the outer tubular housing
d1=4.75*0.0254 # (m) diameter of the inner cold cylinder before
               # cutting any grooves

fig,ax=plt.subplots()

ax.set_xlim([-d2/2,d2/2])
ax.set_ylim([-d2/2,d2/2])

# outer housing
circle1=plt.Circle((0,0),d2/2,color='r',fill=False)
ax.add_artist(circle1)

groove_depth=0.1*0.0254 # m
groove_width=0.06*0.0254 # m
n=ngrooves=124

# shorter names
r=d1/2
w=groove_width
d=groove_depth

for groove in range(ngrooves):
    theta=360./ngrooves
    center_angle=groove*theta

    alpha=center_angle*pi/180
    x=r*cos(alpha)
    y=r*sin(alpha)

    dalpha=asin(w/2/r)

    xsl=r*cos(alpha+dalpha)
    ysl=r*sin(alpha+dalpha)
    xel=xsl-d*cos(alpha)
    yel=ysl-d*sin(alpha)
    
    xsr=r*cos(alpha-dalpha)
    ysr=r*sin(alpha-dalpha)
    xer=xsr-d*cos(alpha)
    yer=ysr-d*sin(alpha)
    
    # draw line at edge of each groove
    line=plt.plot([xsl,xel],[ysl,yel],color='black')
    line=plt.plot([xsr,xer],[ysr,yer],color='black')
    # and bottom of groove
    line=plt.plot([xel,xer],[yel,yer],color='black')

    alpha_deg=alpha*180/pi
    dalpha_deg=dalpha*180/pi
    arc=patches.Arc((0,0),2*r,2*r,0,alpha_deg+dalpha_deg,alpha_deg+theta-dalpha_deg)
    ax.add_patch(arc)

# Calculation of perimeter of all those grooves
pgroove=2*d+w # inner "U" of a groove
parc=(theta-2*dalpha_deg)*pi/180*r # outer arc length between two grooves
perimeter=ngrooves*(pgroove+parc)


# Calculation of area for flow
annulus=pi*(d2**2-d1**2)/4
agroove=d*w # area of one groove # approximately
aeps=((2*dalpha)/(2*pi))*pi*r**2-2*0.5*(w/2)*(r*cos(dalpha)) # area
                                                             # between
                                                             # the arc
                                                             # and the
                                                             # area of
                                                             # the
                                                             # groove
                                                             # on the
                                                             # previous
                                                             # line
agroove_total=agroove+aeps
agrooves=agroove_total*ngrooves
area=annulus+agrooves

plt.show()

a=area # (m^2) area for fluid flow
p=perimeter+pi*d2 # (m) perimeter of flow region
dh=4*a/p # (m) hydraulic diameter
mu=3.5e-5 #Pa*s
kt=0.104 #W/m*K
L=10*0.0254
#
#p_psi150=150. # PSI
#p150=p_psi150*6894.76 # Pa
#
#Tsat150=CP.PropsSI('T','P',p150,'Q',1,fluid) #from cool prop
#
#
#HV150=CP.PropsSI('H','P',p150,'Q',1,fluid)
#
#HL150=CP.PropsSI('H','P',p150,'Q',0,fluid)
#
#ifg150 = HV150 - HL150 #J/kg also from cool prop heat(enthalpy)
#cp=6565. # J/kg*K
#Tw=20.7
#rho=168 #kg/m3
#rho_v150=CP.PropsSI('D','T',Tsat150,'Q',1,fluid) #kg/m3
#ifge150= ifg150 + 0.68*cp*(Tsat150-Tw)
#
#Pr=(mu*cp)/(kt)
#
##jakob number
#
#Ja150 = (cp*(Tsat150-Tw))/ifge150
#
##Rayleigh number
#
#Ra150 = (g*rho*(rho-rho_v150)*L**3*Pr)/(mu**2)
#
#
#LoDcrit150 = 0.007*(Ra150/Ja150)**(1/4)

#print(LoDcrit)
#
#print(L/d1)
#
#print(dh)
#
#Nu150=0.9428*(Ra150/Ja150)**(1/4)
#
##assuming laminar correlation
#
#Re150 = (4*Nu150*Ja150)/Pr
#
#print(Re150)
#
##gratez
#
#Gr150 = (g*rho*(rho-rho_v150)*L**3)/(mu**2)
#
#Nu2150 = ((Pr) / (4*Ja150))*(((Gr150**(1/3)*Ja150) / (0.27*Pr)) + 4.815)**(1/1.22)
#
#Re2150 = (4*Nu2150*Ja150)/Pr
#
#print('Re2 is %f' %Re2150)
#
#hc150 = Nu2150*kt/L

Aw = perimeter*L #(2*n*w+2*n*d)*L (Aw_wet) 3952 compared to 3956
#
#Qdot150 = hc150*Aw*(Tsat150-Tw)
#
#mdot150 = Qdot150/ifge150
#
#
#print('heat coeff is %f W/m2K' %hc150)
#print('convection rate %f W' %Qdot150)
#print('mass flux is %f kg/s'%mdot150 )
#
#print( )
#print()
#print()

p_psi16=16. # PSI
p16=p_psi16*6894.76 # Pa

Tsat16=CP.PropsSI('T','P',p16,'Q',1,fluid) #from cool prop


HV16=CP.PropsSI('H','P',p16,'Q',1,fluid)

HL16=CP.PropsSI('H','P',p16,'Q',0,fluid)

ifg16 = HV16 - HL16 #J/kg also from cool prop heat(enthalpy)
cp=6565. # J/kg*K
Tw=19.8
rho=168 #kg/m3
rho_v16=CP.PropsSI('D','T',Tsat16,'Q',1,fluid) #kg/m3

print('the mass of d2 is %f kg.' %(rho_v16*0.102318))

ifge16= ifg16 + 0.68*cp*(Tsat16-Tw)

Pr=(mu*cp)/(kt)

#jakob number

Ja16 = (cp*(Tsat16-Tw))/ifge16

#Rayleigh number

Ra16 = (g*rho*(rho-rho_v16)*L**3*Pr)/(mu**2)


LoDcrit16 = 0.007*(Ra16/Ja16)**(1/4)

#print(LoDcrit)
#
#print(L/d1)
#
#print(dh)

Nu16=0.9428*(Ra16/Ja16)**(1/4)

#assuming laminar correlation

Re16 = (4*Nu16*Ja16)/Pr

print(Re16)

#gratez

Gr16 = (g*rho*(rho-rho_v16)*L**3)/(mu**2)

Nu216 = ((Pr) / (4*Ja16))*(((Gr16**(1/3)*Ja16) / (0.27*Pr)) + 4.815)**(1/1.22)

Re216 = (4*Nu216*Ja16)/Pr

print('Re2 is %f' %Re216)

hc16 = Nu216*kt/L

Qdot16 = hc16*Aw*(Tsat16-Tw)

mdot16 = Qdot16/ifge16


print('heat coeff is %f W/m2K' %hc16)
print('convection rate %f W' %Qdot16)
print('mass flux is %f kg/s'%mdot16 )
print()


rho_cold=PropsSI('D','T',20,'P',(16*6894.76),fluid)

rho_cold293=PropsSI('D','T',293,'P',(16*6894.76),fluid)

#print(rho_cold)

rho_293=PropsSI('D','T',293,'P',(150*6894.76),fluid)

#print(rho_293)

Vts=0.102318 #m

Vtank = (Vts*(rho_cold-rho_293))*(1/(rho_293-rho_cold293))

print(Vtank)

Tcold=20 #K
Thot=293 #K

Plow=16*6894.76 #Pa
Phigh=150*6894.76 #Pa

Ts = np.linspace(15, 293, 100)
ps = PropsSI('P','T',Ts,'Q',0,'Deuterium')

plt.plot(Ts,ps,'orange')
plt.title('Saturation Curve of Deuterium')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Temp (K)')
plt.text(22.5, 1.25e6, 'liquid', rotation = 45)
plt.text(35, 5e4, 'gas', rotation = 45)
plt.show()

#From there pick Psat off of graph and plug in for p

p= 0.6e6 # Pa

Tsat=CP.PropsSI('T','P',p,'Q',1,fluid) #from cool prop
print()
print('the temp at saturated pressure is %f K' %Tsat)

HV=CP.PropsSI('H','P',p,'Q',1,fluid)

HL=CP.PropsSI('H','P',p,'Q',0,fluid)

ifg = HV - HL #J/kg also from cool prop heat(enthalpy)
cp=6565. # J/kg*K
Tw=19.8
rho=168 #kg/m3
rho_v=CP.PropsSI('D','T',Tsat,'Q',1,fluid) #kg/m3

#print('the mass of d2 is %f kg.' %(rho_v*0.102318))

ifge= ifg + 0.68*cp*(Tsat-Tw)

Pr=(mu*cp)/(kt)

#jakob number

Ja = (cp*(Tsat-Tw))/ifge

#Rayleigh number

Ra = (g*rho*(rho-rho_v)*L**3*Pr)/(mu**2)


LoDcrit = 0.007*(Ra/Ja)**(1/4)

#print(LoDcrit)
#
#print(L/d1)
#
#print(dh)

Nu=0.9428*(Ra/Ja)**(1/4)

#assuming laminar correlation

Re = (4*Nu*Ja)/Pr

print('The Reynolds number is %f' %Re)

#gratez

Gr = (g*rho*(rho-rho_v)*L**3)/(mu**2)

Nu2 = ((Pr) / (4*Ja))*(((Gr**(1/3)*Ja) / (0.27*Pr)) + 4.815)**(1/1.22)

Re2 = (4*Nu2*Ja)/Pr

print('Re2 is %f' %Re2)

hc = Nu2*kt/L

Qdot = hc*Aw*(Tsat-Tw)

mdot = Qdot/ifge


print('heat coeff is %f W/m2K' %hc)
print('convection rate %f W' %Qdot)
print('mass flux is %f kg/s'%mdot)
print()




#
#ps2 = np.linspace(int(Plow), int(Phigh), int(100*6894.76))
#Ts2 = PropsSI('T','P',ps2,'Q',0,'Deuterium')
#
#plt.plot(Ts2,ps2,'orange')
#plt.title('Saturation Curve of Deuterium')
#plt.ylabel('Pressure')
#plt.xlabel('Temp')
#plt.text(22.5, 1.25e6, 'liquid', rotation = 45)
#plt.text(35, 5e4, 'gas', rotation = 45)
#plt.show()
#



