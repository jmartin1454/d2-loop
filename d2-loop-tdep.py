#!/usr/bin/python

# 1D liquid deuterium thermosyphon simulation
# by Jeff and Kiera
# June 19, 2020

from math import *

# gravity
g=9.8 # m/s^2, acceleration due to gravity

# deuterium properties
beta_t=.012 # K^(-1), relative slope of density with temperature (a la
            # Boussinesq)
rho_0=168. # kg/m^3
T_0=21.
def rho(T):
    return rho_0*(1.-beta_t*(T-T_0)) # kg/m^3
cp=6565. # J/(kg-K), specific heat capacity of LD2

# heat input into LD2 moderator vessel
q_mod_total=60 # W, total heat deposited into moderator

# HEX parameters
hc=300. # W/(m^2-K), heat transfer coefficient in HEX
T_cold=19. # K, temperature of the cold wall

# initial conditions
T_initial=T_cold # initial temperature for the whole loop
w=0 # kg/s, mass flow rate


# geometry is
# - a downward-slanted (single spiral) HEX. "hex"
# - a downcomer, which goes straight down "down"
# - a horizontal section "right"
# - a vertical right circular cylinder, with fairly large dimensions, about 125 L "mod"
# - a vertical pipe "rise"
# - another horizontal section going back over to the top of the HEX "left"

n_per=10 # common number of sub-elements per major element

# each major part of the geometry will have a number of constants to
# define it, and a number of arrays containing temperatures and
# geometrical information (lengths along the loop, heights, etc.).  To
# make looping easier, they all get concatenated together after
# setting up.

L_hex=4 # m, length of hex (could be helix)
D_hex=0.015 # m
n_hex=n_per
T_hex=[T_initial]*n_hex
s_hex=[x*L_hex/(n_hex*1.) for x in range(0,n_hex)]
source_hex=[-4*hc*(T_hex[x]-T_cold)/(D_hex*rho_0*cp) for x in range(0,n_hex)]
A_hex=[pi*D_hex**2/4]*n_hex
ds_hex=[L_hex/(n_hex*1.)]*n_hex # step sizes in hex
copper_tube_length=10.*0.0254 # m, physical length of hex
sinalpha=copper_tube_length/L_hex
top_z_hex=2.
z_hex=[top_z_hex-x*L_hex/(n_hex*1.)*sinalpha for x in range(0,n_hex)]

L_down=2 # m, length of downcomer
D_down=0.015 # m, diameter of downcomer
n_down=n_per
T_down=[T_initial]*n_down
s_down=[L_hex+x*L_down/(n_down*1.) for x in range(0,n_down)]
A_down=[pi*D_down**2/4]*n_down
ds_down=[L_down/(n_down*1.)]*n_down
top_z_down=top_z_hex-copper_tube_length
z_down=[top_z_down-x*L_down/(n_down*1.) for x in range(0,n_down)]

L_right=2 # m, length to moderator vessel
D_right=0.015 # m
n_right=n_per
T_right=[T_initial]*n_right
s_right=[L_hex+L_down+x*L_right/(n_right*1.) for x in range(0,n_right)]
A_right=[pi*D_right**2/4]*n_right
ds_right=[L_right/(n_right*1.)]*n_right
bottom_z=top_z_hex-copper_tube_length-L_down
z_right=[bottom_z]*n_right

L_mod=0.5 # m, length of right circular cylinder
D_mod=0.5 # m, diameter
n_mod=n_per
T_mod=[T_initial]*n_mod
s_mod=[L_hex+L_down+L_right+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]
q_h=q_mod_total/(L_mod*pi*D_mod)
source_mod=[4*q_h/(D_mod*rho_0*cp)]*n_mod
A_mod=[pi*D_mod**2/4]*n_mod
ds_mod=[L_mod/(n_mod*1.)]*n_mod
z_mod=[bottom_z+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]

L_rise=top_z_hex-bottom_z-L_mod # m, height of rise
D_rise=0.03 # m, diameter
n_rise=n_per
T_rise=[T_initial]*n_rise
s_rise=[L_hex+L_down+L_right+L_mod+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]
A_rise=[pi*D_rise**2/4]*n_rise
ds_rise=[L_rise/(n_rise*1.)]*n_rise
z_rise=[bottom_z+L_mod+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]

L_left=2 # m, length to moderator vessel
D_left=0.03 # m
n_left=n_per
T_left=[T_initial]*n_left
s_left=[L_hex+L_down+L_right+L_mod+L_rise+x*L_left/(n_left*1.) for x in range(0,n_left)]
A_left=[pi*D_left**2/4]*n_rise
ds_left=[L_left/(n_left*1.)]*n_left
z_left=[top_z_hex]*n_left


# concatenating all the arrays together for easier looping

T_array=T_hex+T_down+T_right+T_mod+T_rise+T_left
n_array=len(T_array)
source_array=source_hex+[0.]*n_down+[0.]*n_right+source_mod+[0.]*n_rise+[0.]*n_left
A_array=A_hex+A_down+A_right+A_mod+A_rise+A_left
s_array=s_hex+s_down+s_right+s_mod+s_rise+s_left
ds_array=ds_hex+ds_down+ds_right+ds_mod+ds_rise+ds_left
z_array=z_hex+z_down+z_right+z_mod+z_rise+z_left
print(n_array,T_array,source_array,A_array,ds_array,z_array)


# sets the beam current, by modifying the appropriate portion of the
# source_array where the moderator is
def set_beam_current(curr):
    power=curr/40.*60. # W, constant of proportionality to power
    q_h=power/(L_mod*pi*D_mod)
    source_mod=[4*q_h/(D_mod*rho_0*cp)]*n_mod
    for i in range(0,n_mod):
        source_array[n_hex+n_down+n_right+i]=source_mod[i]
    return


# main time-stepping loop
n_tsteps=1000000
dt=.1 # s
beam_cycle=240 # s
beam_on=60 # s
for tstep in range(0,n_tsteps):
    t=dt*tstep
    # update temperatures
    for nstep in range(0,n_array):
        dTemp=dt*(-(w/(A_array[nstep]*rho_0))*(T_array[nstep]-T_array[nstep-1])/ds_array[nstep]+source_array[nstep])
        # second derivative term with diffusivity
        # =alpha*(T[i+1]-2T[i]+T[i-1])/(delta s)^2
        T_array[nstep]=T_array[nstep]+dTemp
    # update w
    # rho*dz integral
    rho_integral=0
    for nstep in range(0,n_array):
        rho_integral=rho_integral-rho_0*beta_t*g*T_array[nstep]*(z_array[nstep]-z_array[nstep-1])
    # friction term -- putting this calculation here would allow us
    # eventually to put in temperature/time dependence of friction
    # factor.
    foa2_sum=0.
    f=0.04
    Gamma=0.
    for nstep in range(0,n_array):
        D=(4*A_array[nstep]/pi)**0.5
        # Re=D_H*w/A/mu
        # f=64/Re # for small w, f~1/w -> infty, but w**2*f ~ w -> 0
        # fRe=96 in laminar hex, maybe
        foa2_sum=foa2_sum+f*ds_array[nstep]/D/A_array[nstep]**2
        Gamma=Gamma+ds_array[nstep]/A_array[nstep]
    friction_term=foa2_sum*w**2/(2*rho_0)
    # dw step
    dw=(dt/Gamma)*(-friction_term-rho_integral) # Vijayan (4.25)
    w=w+dw
    if(t%beam_cycle<beam_on):
        set_beam_current(40.)
    else:
        set_beam_current(0.)
    sparse=100 # sparseness of standard output
    if(tstep%sparse==0):
        print('This is time %f and w is %f'%(t,w))
        #print(T_array)
        #print(source_array)
        #print
    for nstep in range(0,n_hex):
        # make hc depend on space and time via correlations
        # Nu=4.8608 or some constant, laminar case
        # hc=Nu*kt/D_h will also be constant
        # in turbulent case, need Re to be calculated.
        source_array[nstep]=-4*hc*(T_array[nstep]-T_cold)/(D_hex*rho_0*cp)
