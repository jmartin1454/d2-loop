#!/usr/bin/python

from math import *

def rho(T):
    return 163.-(T-21.) # kg/m^3
    
    
#rel bt density and temp found by taking the plot of T as a fun of den from coolprop values

T = -2.4506rho + 220.82


T_initial=20.
w=0 # kg/s, mass flow rate
q_mod_total=60 # W, total heat deposited into moderator
cp=6565. # J/(kg-K), specific heat capacity of LD2
hc=300. # W/(m^2-K), heat transfer coefficient in HEX

T_cold=19. # K

L_hex=4 # m, length of hex
D_hex=0.015 # m
n_hex=10
T_hex=[T_initial]*n_hex
s_hex=[x*L_hex/(n_hex*1.) for x in range(0,n_hex)]
source_hex=[-4*hc*(T_hex[x]-T_cold)/(D_hex*rho(21.)*cp) for x in range(0,n_hex)]


L_down=2 # m, length of downcomer
D_down=0.015 # m, diameter of downcomer
n_down=10
T_down=[T_initial]*n_down
s_down=[L_hex+x*L_down/(n_down*1.) for x in range(0,n_down)]

L_right=2 # m, length to moderator vessel
D_right=0.015 # m
n_right=10
T_right=[T_initial]*n_right
s_right=[L_hex+L_down+x*L_right/(n_right*1.) for x in range(0,n_right)]

L_mod=0.5 # m, length of right circular cylinder
D_mod=0.5 # m, diameter
n_mod=10
T_mod=[T_initial]*n_mod
s_mod=[L_hex+L_down+L_right+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]
q_h=q_mod_total/(L_mod*pi*D_mod)
source_mod=[4*q_h/(D_mod*rho(21.)*cp)]*n_mod

L_rise=1.6 # m, height of rise
D_rise=0.03 # m, diameter
n_rise=10
T_rise=[T_initial]*n_rise
s_rise=[L_hex+L_down+L_right+L_mod+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]

L_left=2 # m, length to moderator vessel
D_left=0.03 # m
n_left=10
T_left=[T_initial]*n_left
s_left=[L_hex+L_down+L_right+L_mod+L_rise+x*L_left/(n_left*1.) for x in range(0,n_left)]


T_array=T_hex+T_down+T_right+T_mod+T_rise+T_left
n_array=len(T_array)
source_array=source_hex+[0.]*n_down+[0.]*n_right+source_mod+[0.]*n_rise+[0.]*n_left
print(n_array,T_array,source_array)

n_tsteps=100
dt=1.0 # s
for tstep in range(0,n_tsteps):
    t=dt*tstep




















