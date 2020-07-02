#!/usr/bin/python

from math import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.misc import derivative


kt=0.104

mu=3.5e-5

g=9.8 # m/s^2
beta_t=.012 # K^(-1), relative slope of density with temperature (a la
            # Boussinesq)
rho_0=168. # kg/m^3
T_0=21.
def rho(T):

    return 220.82-(T/2.4506) # kg/m^3
    
    
#rel bt density and temp found by taking the plot of T as a fun of den from coolprop values

#T = -2.4506rho + 220.82

#    return rho_0*(1.-beta_t*(T-T_0)) # kg/m^3

#k values

#first sudden cont
A1= (pi*0.127**2)/4 #m^2 area of pipe bf exp  L = 0.18641 D = 0.127

A3= (pi*0.03808**2)/4 #m^2 area of pipe after exp  L = 0.02093 D = 0.03808

Kcont=(1/2)*(1-A3/A1)**(3/4)

# For second sudden cont

A12= (pi*0.03808**2)/4 #m^2 area of pipe bf exp  L = 0.02093 D = 0.03808

A32= (pi*0.0127**2)/4 #m^2 area of pipe after exp L =  D = 0.0127

Kcont2=(1/2)*(1-A3/A1)**(3/4)

K45 = 0.3

# For sudden expansion from hydrualic_Resistance.pdf for turb??

Aexp1= (pi*0.0127**2)/4 #m^2 area of pipe bf exp

Aexp2= (pi*((159.5+29.5)*0.0254)**2) #m^2 area of pipe after exp

Kexp=(1-Aexp1/Aexp2)**2


K180=2.2


# For sudden contraction from hydrualic_Resistance.pdf

A1=Aexp2 #m^2 before cont

Acont3= (pi*0.03175**2)/4 #m^2 after

Kcont3=(1/2)*(1-A3/A1)**(3/4)

#90 deg

D= 0.03175 # m
A=pi*D**2/4 # m^2

K902=0.9

#valve ? sudden exp ?

A1= (pi*0.03175**2)/4 #m^2 area of pipe bf exp

A2= (pi*0.127**2)/4 #m^2 area of pipe after exp

Kvalve=10

# For sudden expansion from hydrualic_Resistance.pdf for turb??

A1= (pi*0.03175**2)/4 #m^2 area of pipe bf exp

A2= (pi*0.127**2)/4 #m^2 area of pipe after exp

Kexp2=(1-A1/A2)**2


w=0 # kg/s, mass flow rate
q_mod_total=60 # W, total heat deposited into moderator
cp=6565. # J/(kg-K), specific heat capacity of LD2
hc=300. # W/(m^2-K), heat transfer coefficient in HEX

T_cold=19.8 # K
T_initial=T_cold

n_per=10

L_hex=4 # m, length of hex (could be helix)
D_hex=0.015 # m
n_hex=n_per
T_hex=[T_initial]*n_hex
s_hex=[x*L_hex/(n_hex*1.) for x in range(0,n_hex)]
source_hex=[-4*hc*(T_hex[x]-T_cold)/(D_hex*rho(21.)*cp) for x in range(0,n_hex)]
A_hex=[pi*D_hex**2/4]*n_hex
ds_hex=[L_hex/(n_hex*1.)]*n_hex # step sizes in hex
copper_tube_length=10.*0.0254 # m, physical length of hex
sinalpha=copper_tube_length/L_hex
top_z_hex=2.
z_hex=[top_z_hex-x*L_hex/(n_hex*1.)*sinalpha for x in range(0,n_hex)]
P_hex=[pi*D_hex]*n_hex



L_down=1.937 # m, length of downcomer
D_down=0.0134 # m, diameter of downcomer
n_down=n_per
T_down=[T_initial]*n_down
s_down=[L_hex+x*L_down/(n_down*1.) for x in range(0,n_down)]
A_down=[pi*D_down**2/4]*n_down
A2_down=pi*D_down**2/4
ds_down=[L_down/(n_down*1.)]*n_down
top_z_down=top_z_hex-copper_tube_length
z_down=[top_z_down-x*L_down/(n_down*1.) for x in range(0,n_down)]
P_down=[pi*D_down]*n_down



L_right=(1.27+0.8) # m, length to moderator vessel
D_right=0.0134 # m
n_right=n_per
T_right=[T_initial]*n_right
s_right=[L_hex+L_down+x*L_right/(n_right*1.) for x in range(0,n_right)]
A_right=[pi*D_right**2/4]*n_right
ds_right=[L_right/(n_right*1.)]*n_right
bottom_z=top_z_hex-copper_tube_length-L_down
z_right=[bottom_z]*n_right
P_right=[pi*D_right]*n_right



L_mod=0.5 # m, length of right circular cylinder
D_mod=0.5 # m, diameter
n_mod=n_per
T_mod=[T_initial]*n_mod
s_mod=[L_hex+L_down+L_right+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]
q_h=q_mod_total/(L_mod*pi*D_mod)
source_mod=[4*q_h/(D_mod*rho(21.)*cp)]*n_mod
A_mod=[pi*D_mod**2/4]*n_mod
ds_mod=[L_mod/(n_mod*1.)]*n_mod
z_mod=[bottom_z+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]
P_mod=[pi*D_mod]*n_mod



L_rise=1.5 #top_z_hex-bottom_z-L_mod # m, height of rise
D_rise=0.03175 # m, diameter
n_rise=n_per
T_rise=[T_initial]*n_rise
s_rise=[L_hex+L_down+L_right+L_mod+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]
A_rise=[pi*D_rise**2/4]*n_rise
ds_rise=[L_rise/(n_rise*1.)]*n_rise
z_rise=[bottom_z+L_mod+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]
P_rise=[pi*D_rise]*n_rise



L_left=1.75 # m, length to moderator vessel
D_left=0.03175 # m
n_left=n_per
T_left=[T_initial]*n_left
s_left=[L_hex+L_down+L_right+L_mod+L_rise+x*L_left/(n_left*1.) for x in range(0,n_left)]
A_left=[pi*D_left**2/4]*n_left
ds_left=[L_left/(n_left*1.)]*n_left
z_left=[top_z_hex]*n_left
P_left=[pi*D_left]*n_left



T_array=T_hex+T_down+T_right+T_mod+T_rise+T_left
n_array=len(T_array)
source_array=source_hex+[0.]*n_down+[0.]*n_right+source_mod+[0.]*n_rise+[0.]*n_left
A_array=A_hex+A_down+A_right+A_mod+A_rise+A_left
P_array=P_hex+P_down+P_right+P_rise+P_mod+P_left
s_array=s_hex+s_down+s_right+s_mod+s_rise+s_left
ds_array=ds_hex+ds_down+ds_right+ds_mod+ds_rise+ds_left
ds2_array=ds_hex+ds_down+ds_right+ds_rise+ds_left
z_array=z_hex+z_down+z_right+z_mod+z_rise+z_left
#print(n_array,T_array,source_array,A_array,ds_array,z_array)



wvalue=[]

tvalue=[]

Tminvalue=[]

Tmaxvalue=[]

Tend_array=[]

T1=[]

fvalue=[]

Revalue=[]

# sets the beam current, by modifying the appropriate portion of the
# source_array where the moderator is
def set_beam_current(curr):
    power=curr/40.*60. # W, constant of proportionality to power
    q_h=power/(L_mod*pi*D_mod)
    source_mod=[4*q_h/(D_mod*rho_0*cp)]*n_mod
    for i in range(0,n_mod):
        source_array[n_hex+n_down+n_right+i]=source_mod[i]
    return


Nu=4.8608 #or some constant, laminar case
hc=Nu*kt/D #will also be constant
n_tsteps=6000000
dt=.0001 # s
# consider adaptive time steps see Vijayan eq. (4.99)
beam_cycle=240 # s
beam_on=60 # s
alpha=kt/(rho_0*cp) #J/smk kgm2/ss2mK    kgm/s3K * 1/kg/m3 * 1/J/kgK
# kgm m3 kg K / s3 K kg J     kg m4 s2/ s3 kg m2 m2 /s
for tstep in range(0,n_tsteps):
    t=dt*tstep
    temp = interpolate.interp1d(ds_array, T_array)
    for nstep in range(0,n_array):
        dTemp=dt*(-(w/(A_array[nstep]*rho(21.)))*(T_array[nstep]-T_array[nstep-1])/ds_array[nstep]+source_array[nstep] + alpha*(T_array[nstep]-2*T_array[nstep-1]+T_array[nstep-2])/(ds_array[nstep]**2))
        T_array[nstep]=T_array[nstep]+dTemp
#        def Temperature(x):
#            return temp(x)
#        def Temperature_prime(x):
#            return derivative(Temperature, x, 1e-6)
#        def Temperature_primeprime(x):
#            return derivative(Temperature_prime, x, 1e-6)
#        dTemp=dt*(-(w/(A_array[nstep]*rho(21.)))*Temperature_prime + source_array[nstep] + alpha*Temperature_primeprime ) #alpha*(T_array[nstep]-2*T_array[nstep-1]+T_array[nstep-2]/(ds_array**2))
#        T_array[nstep]=T_array[nstep]+dTemp
        
    # update temperatures
#    for nstep in range(0,n_array):
    # update w
    # rho integral
    rho_integral=0
    for nstep in range(0,n_array):
        rho_integral=rho_integral-rho_0*beta_t*g*T_array[nstep]*(z_array[nstep]-z_array[nstep-1])
    # friction term
    foa2_sum=0.
    Gamma=0.
    #f=0.03 #f_hex+f_down+f_left+f_right+f_rise
    #sumfLoDA2=(((f_hex*L_hex)/(D_hex*((pi*(D_hex/2)**2))**2)))+(((f_down*L_down)/(D_down*((pi*(D_down/2)**2))**2)))+(((f_right*L_right)/(D_right*((pi*(D_right/2)**2))**2)))+(((f_rise*L_rise + 2*K45)/(D_rise*((pi*(D_rise/2)**2))**2)))+(((f_left*L_left)/(D_left*((pi*(D_left/2)**2))**2))) + (Kcont )*(1/A3**2) +(Kcont2 )*(1/A32**2)+(Kexp)*(1/Aexp2**2)+(Kcont3)*(1/Acont3**2)+(K902 )*(1/A**2)+(Kvalve )*(1/A2**2)+(Kexp2 )*(1/A2**2)
    for nstep in range(0,n_array):
        D=(4*A_array[nstep]/pi)**0.5
        D_h=4*A_array[nstep]/P_array[nstep]
        Re=D_h*w/(A_array[nstep]*mu)
        Revalue.append(Re)
        #f=64/Re # for small w, f~1/w -> infty, but w**2*f ~ w -> 0
        # fRe=96 in laminar hex, maybe
        if w < 0.0000003:
            f = 0.03
        else :
            f=64/Re
        #print('fRe=%f'%(f*Re))
        fvalue.append(f)
        foa2_sum=foa2_sum + f*ds_array[nstep]/(D_h*A_array[nstep]**2 )
        Gamma=Gamma+ds_array[nstep]/A_array[nstep]
    friction_term=foa2_sum*w**2/(2*rho_0)
    # dw step
    dw=(dt/Gamma)*(-friction_term-rho_integral) # Vijayan (4.25)
    w=w+dw
    if(t%beam_cycle<beam_on):
        set_beam_current(10.)
    else:

        set_beam_current(10.)
    sparse=100000 # sparseness of standard output
    if(tstep%sparse==0):
        print('This is time %f and w is %f'%(t,w))
        print(min(T_array),max(T_array))
        #print(T_array)
        Tminvalue.append(min(T_array))
        Tmaxvalue.append(max(T_array))
        wvalue.append(w)
        tvalue.append(t)
        #print(source_array)
        #print
    for nstep in range(0,n_hex):
        source_array[nstep]=-4*hc*(T_array[nstep]-T_cold)/(D_hex*rho(21.)*cp)



plt.plot(tvalue,wvalue,'r:')
plt.ylabel('w (kg/s)')
plt.title('Mass Flux as a Function of Time')
plt.xlabel('Time (s)')
plt.show()



#plt.plot(Revalue,fvalue,'ro')
#plt.xlabel('Reynolds number')
#plt.ylabel('fricition factor')
#plt.yscale('log')
#plt.xscale('log')
#plt.show()




plt.plot(tvalue,Tminvalue,'b:',label='Min Temp')
plt.plot(tvalue,Tmaxvalue,'r:',label='Max Temp')
plt.legend(loc='upper left')
plt.title('Temperature as a Function of Time')
plt.ylabel('Temperature (K)')
plt.xlabel('Time (s)')
plt.show()



plt.plot(T_array,'b:')
#plt.plot(T1,'r:')
plt.ylabel('Temperature (K)')
plt.xlabel('Position')
plt.show()


plt.plot(s_array,T_array,'b:')
#plt.plot(T1,'r:')
plt.title('Temperature as a Function of Physical Position Around the Loop')
plt.ylabel('Temperature (K)')
plt.xlabel('Position (m)')
plt.show()

plt.plot(z_array,T_array,'b:')
#plt.plot(T1,'r:')
plt.ylabel('Temperature (K)')
plt.xlabel('z_array (m)')
plt.show()

print(foa2_sum)
#print(T1)


# make hc depend on space and time via correlations
Nu=4.8608 #or some constant, laminar case
hc=Nu*kt/D #will also be constant
# in turbulent case, need Re to be calculated.
source_array[nstep]=-4*hc*(T_array[nstep]-T_cold)/(D_hex*rho_0*cp)


























