#!/usr/bin/python3

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy import interpolate
from scipy.misc import derivative
from numpy import asarray
from numpy import savetxt

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-c", "--camhex", dest="camhex",
                  default=False, action="store_true",
                  help="Use original Cam HEX design")

(options, args) = parser.parse_args()


kt=0.104
Nu=4.8608
mu=3.5e-5

##############################################
# material parameters for some walls
# https://trc.nist.gov/cryogenics/calculators/propcalc.html
##############################################
# 304 stainless steel
rho_steel=7850 # kg/m^3
kt_steel=2.169 # W/(m*K)
cp_steel=13.452 # J/(kg*K)
# 6061-T6 Al
rho_al=2700 # kg/m^3
kt_al=28.428 # W/(m*K)
cp_al=8.854 # J/(kg*K)


g=9.8 # m/s^2
beta_t=.012 # K^(-1), relative slope of density with temperature (a la
            # Boussinesq)
rho_0=168. # kg/m^3
T_0=21.
def rho(T):

#    return 220.82-(T/2.4506) # kg/m^3

#rel bt density and temp found by taking the plot of T as a fun of den from coolprop values

#T = -2.4506rho + 220.82

    return rho_0*(1.-beta_t*(T-T_0)) # kg/m^3

# Note: using any version of rho(T) other than the line above is
# inconsistent with equations below which also rely on beta_t.  If you
# want to use any other version, you must update beta_t and rho_0
# rather than hard-coding your own Boussinesq formula.

#k values

#first sudden cont
A1= (pi*0.127**2)/4 #m^2 area of pipe bf exp  L = 0.18641 D = 0.127

A3= (pi*0.03808**2)/4 #m^2 area of pipe after exp  L = 0.02093 D = 0.03808

Kcont=(1/2)*(1-A3/A1)**(3/4)

KoA2cont1=Kcont/(A3**2)

# For second sudden cont

A12= (pi*0.03808**2)/4 #m^2 area of pipe bf exp  L = 0.02093 D = 0.03808

A32= (pi*0.0127**2)/4 #m^2 area of pipe after exp L =  D = 0.0127

Kcont2=(1/2)*(1-A3/A1)**(3/4)

KoA2cont2=Kcont2/A32**2

# 45 deg turn

K45 = 0.3

D_down=0.0134

KoA245=K45/(pi*D_down**2/4)**2

# For sudden expansion from hydrualic_Resistance.pdf for turb??

Aexp1= (pi*0.0127**2)/4 #m^2 area of pipe bf exp

Aexp2= (pi*((159.5+29.5)*0.0254)**2) #m^2 area of pipe after exp

Kexp=(1-Aexp1/Aexp2)**2

KoA2exp1=Kexp/Aexp2**2

#180 deg turn


K180=2.2

D_mod=0.5

KoA2180=K180/(pi*D_mod**2/4)**2


# For sudden contraction from hydrualic_Resistance.pdf

A1=Aexp2 #m^2 before cont

Acont3= (pi*0.03175**2)/4 #m^2 after

Kcont3=(1/2)*(1-A3/A1)**(3/4)

KoA2cont3=Kcont3/Acont3**2

#90 deg

D= 0.03175 # m
A90=pi*D**2/4 # m^2

K902=0.9

KoA290=K902/A90**2

#valve ? sudden exp ?

A1valve= (pi*0.03175**2)/4 #m^2 area of pipe bf exp

A2valve= (pi*0.127**2)/4 #m^2 area of pipe after exp

Kvalve=10

KoA2valve=Kvalve/A2valve**2

# For sudden expansion from hydrualic_Resistance.pdf for turb??

A1exp2= (pi*0.03175**2)/4 #m^2 area of pipe bf exp

A2exp2= (pi*0.127**2)/4 #m^2 area of pipe after exp

Kexp2=(1-A1exp2/A2exp2)**2

KoA2exp2=Kexp2/A2exp2**2


w=0 # kg/s, mass flow rate
q_mod_total=60/4 # W, total heat deposited into moderator
cp=6565. # J/(kg-K), specific heat capacity of LD2
 # W/(m^2-K), heat transfer coefficient in HEX

T_cold=19.8 # K
T_initial=T_cold+0.1

n_per=10

##############################################
# hex geometry of inital laminar
##############################################

if options.camhex:

    L_hex=Lprime=10*0.0254 # m, length of hex (could be helix)
    d2=4.76*0.0254 # (m) inner diameter of the outer tubular housing
    d1=4.75*0.0254 # (m) diameter of the inner cold cylinder before
    groove_depth=0.1*0.0254 # m
    groove_width=0.06*0.0254 # m
    n=ngrooves=124
    fig,ax=plt.subplots()
    ax.set_xlim([-d2/2,d2/2])
    ax.set_ylim([-d2/2,d2/2])
    circle1=plt.Circle((0,0),d2/2,color='r',fill=False)
    ax.add_artist(circle1)
    r=d1/2
    wd=groove_width
    d=groove_depth
    theta=360./ngrooves
    for groove in range(ngrooves):
        center_angle=groove*theta

        alpha=center_angle*pi/180
        x=r*cos(alpha)
        y=r*sin(alpha)

        dalpha=asin(wd/2/r)

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
    plt.show()
    pgroove=(2*d)+wd # inner "U" of a groove
    parc=(theta-2*dalpha_deg)*pi/180*r # outer arc length between two grooves
    Phc=perimeter=ngrooves*(pgroove+parc)
    annulus=pi*(d2**2-d1**2)/4
    agroove=d*wd # area of one groove # approximately
    aeps=((2*dalpha)/(2*pi))*pi*r**2-2*0.5*(wd/2)*(r*cos(dalpha))
    agroove_total=agroove+aeps
    agrooves=agroove_total*ngrooves
    #D_hex=0.001523 # m hydraulic diameter
    n_hex=n_per
    A=annulus+agrooves
    P=perimeter+(pi*d2)
    P_hex=[perimeter+pi*d2]*n_hex # flow perimeter
    A_hex=[annulus+agrooves]*n_hex
    D_hex=4*A/P
    T_hex=[T_initial]*n_hex
    hc=Nu*kt/D_hex #331.916164 # W/m^2K
    s_hex=[x*L_hex/(n_hex*1.) for x in range(0,n_hex)]
    source_hex=[-4*hc*(T_hex[x]-T_cold)/(D_hex*rho_0*cp) for x in range(0,n_hex)]
    ds_hex=[L_hex/(n_hex*1.)]*n_hex # step sizes in hex
    top_z_hex=(1.937+(10*0.0254))
    z_hex=[top_z_hex-x*L_hex/(n_hex*1.) for x in range(0,n_hex)]
    #Grhex= (g*beta_t*rho**2*Dh**3*(q_h*(L_hex+L_down)/(A*mu*cp)))/(mu**2)
    x_hex=[0.0]*n_hex

    ##############################################
    # HEX stainless wall parameters
    ##############################################
    T_hex_wall=[T_initial]*n_hex
    alpha_hex_wall=[kt_steel/(cp_steel*rho_steel)]*n_hex
    rhocp_hex_wall=[cp_steel*rho_steel]*n_hex
    ID_hex_wall=[d2]*n_hex
    thick_hex_wall=[(5*0.0254-d2)/2]*n_hex # wall thickness, guess 5" OD

else:
    ##############################################
    # hex geometry of helix, taken from https://github.com/jmartin1454/d2-hex/blob/master/Helix-as-tube-v2.ipynb
    ##############################################

    L=L_hex=10*0.0254 #m physical height of tube, from CM.LD2.015.R2
    Ngrooves=1 # number of grooves, the optimal groove will be explored in another notebook 
    D=4.76*0.0254 #m diameter of tube from CM.LD2.015.R2
    wprime= 0.015 #m width of groove
    uprime= 0.01 # m width between grooves
    depth=0.015 # m depth of groove
    sinalpha=(Ngrooves*(wprime+uprime))/(pi*D) #pitch angle
    alpha=asin(sinalpha)

    Lprime=L/sinalpha #m, length of wound groove
    appturns=Lprime/(pi*D) # note:  a small-angle approximation has been used here
    turns=L/(tan(alpha)*pi*D) # correct calculation with no small-angle approximation
    ahelix=Ngrooves*wprime*depth # m^2, total flow area of all the grooves
    phelix=Ngrooves*(2*depth+2*wprime) # m, total flow perimeter of all the grooves 
    Dh=4*ahelix/phelix #m
    pcold=perimeter=Ngrooves*(wprime+2*depth) # m, cold perimeter

    n_hex=n_per
    A=ahelix
    P=phelix
    P_hex=[phelix]*n_hex # flow perimeter
    A_hex=[ahelix]*n_hex
    D_hex=4*A/P
    T_hex=[T_initial]*n_hex
    hc=Nu*kt/D_hex #331.916164 # W/m^2K
    s_hex=[x*Lprime/(n_hex*1.) for x in range(0,n_hex)]
    source_hex=[-4*hc*(T_hex[x]-T_cold)/(D_hex*rho_0*cp)*pcold/phelix for x in range(0,n_hex)]
    ds_hex=[Lprime/(n_hex*1.)]*n_hex # step sizes in hex
    top_z_hex=(1.937+(10*0.0254))
    z_hex=[top_z_hex-x*L/(n_hex*1.) for x in range(0,n_hex)]
    #Grhex= (g*beta_t*rho**2*Dh**3*(q_h*(L_hex+L_down)/(A*mu*cp)))/(mu**2)
    x_hex=[0.0]*n_hex

    ##############################################
    # HEX stainless wall parameters
    ##############################################
    T_hex_wall=[T_initial]*n_hex
    alpha_hex_wall=[kt_steel/(cp_steel*rho_steel)]*n_hex
    rhocp_hex_wall=[cp_steel*rho_steel]*n_hex
    ID_hex_wall=[D]*n_hex
    thick_hex_wall=[(5*0.0254-D)/2]*n_hex # wall thickness, guess 5" OD

    
##############################################
# 1st downcomer geometry
##############################################
L_down=1.937 # m, length of downcomer
D_down=0.0134 # m, diameter of downcomer
n_down=n_per
T_down=[T_initial]*n_down
s_down=[Lprime+x*L_down/(n_down*1.) for x in range(0,n_down)]
Adown=pi*D_down**2/4
A_down=[pi*D_down**2/4]*n_down
ds_down=[L_down/(n_down*1.)]*n_down
top_z_down=top_z_hex-L_hex
z_down=[top_z_down-x*L_down/(n_down*1.) for x in range(0,n_down)]
P_down=[pi*D_down]*n_down
#Grdown = (g*beta_t*rho**2*D_down**3*(q_h*(L_hex+L_down)/(Adown*mu*cp)))/(mu**2)
x_down=[0.0]*n_down

##############################################
# 1st downcomer wall parameters
##############################################
T_down_wall=[T_initial]*n_down
alpha_down_wall=[kt_steel/(cp_steel*rho_steel)]*n_down
rhocp_down_wall=[cp_steel*rho_steel]*n_down
ID_down_wall=[D_down]*n_down
# guess schedule 10, with wall thickness of 0.065 inches
thick_down_wall=[0.065*0.0254]*n_down # m

##############################################
# right section going over to moderator
##############################################
L_right=(1.27+0.8) # m, length to moderator vessel
D_right=0.0134 # m
Aright=pi*D_right**2/4
n_right=n_per
T_right=[T_initial]*n_right
s_right=[Lprime+L_down+x*L_right/(n_right*1.) for x in range(0,n_right)]
A_right=[pi*D_right**2/4]*n_right
ds_right=[L_right/(n_right*1.)]*n_right
bottom_z=top_z_hex-L_hex-L_down
z_right=[bottom_z]*n_right
P_right=[pi*D_right]*n_right
#Grright = (g*beta_t*rho**2*D_right**3*(q_h*(L_hex+L_down)/(Aright*mu*cp)))/(mu**2)
x_right=[x*L_right/(n_right*1.) for x in range(0,n_right)]

##############################################
# right going wall parameters
##############################################
T_right_wall=[T_initial]*n_right
alpha_right_wall=[kt_steel/(cp_steel*rho_steel)]*n_right
rhocp_right_wall=[cp_steel*rho_steel]*n_right
ID_right_wall=[D_right]*n_right
# guess schedule 10, with wall thickness of 0.065 inches
thick_right_wall=[0.065*0.0254]*n_right # m

##############################################
# 2nd downcomer within moderator volume
##############################################
L_down2=0.5 # m, length to moderator vessel
D_down2=D_right # m
Adown2=pi*D_down2**2/4
n_down2=n_per
T_down2=[T_initial]*n_down2
s_down2=[Lprime+L_down+L_right+x*L_down2/(n_down2*1.) for x in range(0,n_down2)]
A_down2=[pi*D_down2**2/4]*n_down2
ds_down2=[L_down2/(n_down2*1.)]*n_down2
bottom_z=top_z_hex-L_hex-L_down
z_down2=[bottom_z-x*L_down2/(n_down2*1.) for x in range(0,n_down2)]
very_bottom_z=bottom_z-L_down2
P_down2=[pi*D_down2]*n_down2
#Grdown2 = (g*beta_t*rho**2*D_down2**3*(q_h*(L_hex+L_down)/(Adown2*mu*cp)))/(mu**2)
x_down2=[L_right]*n_down2


##############################################
# 2nd downcomer wall parameters
##############################################
T_down2_wall=[T_initial]*n_down2
alpha_down2_wall=[kt_al/(cp_al*rho_al)]*n_down2
rhocp_down2_wall=[cp_al*rho_al]*n_down2
ID_down2_wall=[D_down2]*n_down2
# guess schedule 10, with wall thickness of 0.065 inches
thick_down2_wall=[0.065*0.0254]*n_down2 # m

##############################################
# moderator volume
##############################################
L_mod=0.5 # m, length of right circular cylinder
D_mod=0.5 # m, diameter
n_mod=n_per
T_mod=[T_initial]*n_mod
s_mod=[Lprime+L_down+L_right+L_down2+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]
q_h=q_mod_total/(L_mod*pi*D_mod) # W/m^2
source_mod=[4*q_h/(D_mod*rho_0*cp)]*n_mod # W/(J/K)=K/s
A_mod=[pi*D_mod**2/4]*n_mod
Amod=pi*D_mod**2/4
ds_mod=[L_mod/(n_mod*1.)]*n_mod
z_mod=[very_bottom_z+x*L_mod/(n_mod*1.) for x in range(0,n_mod)]
P_mod=[pi*D_mod]*n_mod
#Grmod = (g*beta_t*rho**2*D_mod**3*(q_h*(L_hex+L_down)/(Amod*mu*cp)))/(mu**2)
x_mod=[L_right]*n_mod
hc_mod = Nu*kt/D_mod

##############################################
# moderator wall parameters
##############################################
T_mod_wall=[T_initial]*n_mod
alpha_mod_wall=[kt_al/(cp_al*rho_al)]*n_mod
rhocp_mod_wall=[cp_al*rho_al]*n_mod
ID_mod_wall=[D_mod]*n_mod
# guess schedule 10, with wall thickness of 0.065 inches
thick_mod_wall=[0.065*0.0254]*n_mod # m

##############################################
# calculate total heat deposited to down2 based on proportional heat
# in moderator
##############################################
# make it proportional to moderator total heat, by volume
Q_down2=q_mod_total*(pi*D_down2**2*L_down2/4)/(pi*D_mod**2*L_mod/4) # W
# heat per unit volume in either the moderator or this down2 tube
Q_per_volume=q_mod_total/(pi*D_mod**2*L_mod/4)

source_mod_down2=[Q_per_volume/(rho_0*cp)]*n_down2
#source_mod_down2=[0.]*n_down2
# (W/m^3)/((kg/m^3)*(J/(kg*K)))
# K/s
# used to artificially reduce heat in down2 section
#source_mod_down2[0]=0
#source_mod_down2[1]=0
#source_mod_down2[2]=0

#print(Q_down2)
#print(Q_down2/(pi*D_down2**2*L_down2/4))
#print(Q_per_volume)
#print(4*q_h/D_mod)


##############################################
# left section leaving moderator
##############################################

L_left=L_right # m, length to moderator vessel
D_left=0.03175 # m
n_left=n_per
T_left=[T_initial]*n_left
s_left=[Lprime+L_down+L_right+L_down2+L_mod+x*L_left/(n_left*1.) for x in range(0,n_left)]
A_left=[pi*D_left**2/4]*n_left
Aleft=pi*D_left**2/4
ds_left=[L_left/(n_left*1.)]*n_left
z_left=[bottom_z]*n_left
P_left=[pi*D_left]*n_left
#Grleft = (g*beta_t*rho**2*D_left**3*(q_h*(L_hex+L_down)/(Aleft*mu*cp)))/(mu**2)
x_left=[L_right-x*L_left/(n_left*1.) for x in range(0,n_left)]

##############################################
# left wall
##############################################
T_left_wall=[T_initial]*n_left
alpha_left_wall=[kt_steel/(cp_steel*rho_steel)]*n_left
rhocp_left_wall=[cp_steel*rho_steel]*n_left
ID_left_wall=[D_left]*n_left
# guess schedule 10, with wall thickness of 0.065 inches
thick_left_wall=[0.065*0.0254]*n_left # m

##############################################
# riser
##############################################
L_rise=L_hex+L_down #top_z_hex-bottom_z-L_mod # m, height of rise
D_rise=0.03175 # m, diameter
n_rise=n_per
T_rise=[T_initial]*n_rise
s_rise=[Lprime+L_down+L_right+L_down2+L_mod+L_left+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]
A_rise=[pi*D_rise**2/4]*n_rise
Arise=pi*D_rise**2/4
ds_rise=[L_rise/(n_rise*1.)]*n_rise
z_rise=[bottom_z+x*L_rise/(n_rise*1.) for x in range(0,n_rise)]
P_rise=[pi*D_rise]*n_rise
#Grrise = (g*beta_t*rho**2*D_rise**3*(q_h*(L_hex+L_down)/(Arise*mu*cp)))/(mu**2)
x_rise=[L_right-L_left]*n_rise

##############################################
# riser wall
##############################################
T_rise_wall=[T_initial]*n_rise
alpha_rise_wall=[kt_steel/(cp_steel*rho_steel)]*n_rise
rhocp_rise_wall=[cp_steel*rho_steel]*n_rise
ID_rise_wall=[D_rise]*n_rise
# guess schedule 10, with wall thickness of 0.065 inches
thick_rise_wall=[0.065*0.0254]*n_rise # m

##############################################
# Concatenate arrays to make more easily iterable
##############################################
T_array=np.array(T_hex+T_down+T_right+T_down2+T_mod+T_left+T_rise)
n_array=len(T_array)
A_array=np.array(A_hex+A_down+A_right+A_down2+A_mod+A_left+A_rise)
P_array=np.array(P_hex+P_down+P_right+P_down2+P_mod+P_left+P_rise)
s_array=np.array(s_hex+s_down+s_right+s_down2+s_mod+s_left+s_rise)
z_array=np.array(z_hex+z_down+z_right+z_down2+z_mod+z_left+z_rise)
ds_array=np.array(ds_hex+ds_down+ds_right+ds_down2+ds_mod+ds_left+ds_rise)
x_array=np.array(x_hex+x_down+x_right+x_down2+x_mod+x_left+x_rise)

##############################################
# Concatenate source arrays
##############################################
source_array=source_hex+[0.]*n_down+[0.]*n_right+source_mod_down2+source_mod+[0.]*n_left+[0.]*n_rise
source_fluid_to_wall=[0.]*(n_hex+n_down+n_right+n_down2+n_mod+n_left+n_rise)
source_wall_to_fluid=[0.]*(n_hex+n_down+n_right+n_down2+n_mod+n_left+n_rise)

##############################################
# Concatenate wall arrays
##############################################
T_wall=np.array(T_hex_wall+T_down_wall+T_right_wall+T_down2_wall+T_mod_wall+T_left_wall+T_rise_wall)
alpha_wall=np.array(alpha_hex_wall+alpha_down_wall+alpha_right_wall+alpha_down2_wall+alpha_mod_wall+alpha_left_wall+alpha_rise_wall)
rhocp_wall=np.array(rhocp_hex_wall+rhocp_down_wall+rhocp_right_wall+rhocp_down2_wall+rhocp_mod_wall+rhocp_left_wall+rhocp_rise_wall)
ID_wall=np.array(ID_hex_wall+ID_down_wall+ID_right_wall+ID_down2_wall+ID_mod_wall+ID_left_wall+ID_rise_wall)
thick_wall=np.array(thick_hex_wall+thick_down_wall+thick_right_wall+thick_down2_wall+thick_mod_wall+thick_left_wall+thick_rise_wall)

#Lt=L_hex+L_down+L_right+L_mod+L_left+L_rise
#Dr=(1/Lt)*(Dh*L_hex + D_down*L_down + D_right*L_right + D_mod*L_mod + D_left*L_left + D_rise*L_rise)
#Ng=(Lt/Dr)*()


Tot_vol = (A*Lprime) + (pi*D_down**2*L_down/4) + (pi*D_right**2*L_right/4) + (pi*D_down2**2*L_down2/4) + (pi*D_mod**2*L_mod/4) + (pi*D_left**2*L_left/4) + (pi*D_rise**2*L_rise/4)
print('total vol is %f m3' %Tot_vol )


# graphs of geometry
plt.title('The Cross Section of the Heat Exchanger.')
plt.plot(s_array,z_array,'b:')
plt.title('Loop')
plt.ylabel('Position (m)')
plt.xlabel('Position (m)')
plt.show()

plt.plot(x_array,z_array,'b:')
plt.title('Loop')
plt.ylabel('Vertical Position (m)')
plt.xlabel('Horizontal Position (m)')
plt.show()


wvalue=[]
tvalue=[]
Tminvalue=[]
Tmaxvalue=[]
Tend_array=[]
T1=[]
fvalue=[]


# sets the beam current, by modifying the appropriate portion of the
# source_array where the moderator is
def set_beam_current(curr):
    fractional_power=curr/10. # dimensionless, normalized to 10 uA
    source=Q_per_volume*fractional_power/(rho_0*cp)
    # heat in 2nd downcomer
    source_array[n_hex+n_down+n_right:n_hex+n_down+n_right+n_down2]=[source]*n_down2
    # heat in moderator volume itself
    source_array[n_hex+n_down+n_right+n_down2:n_hex+n_down+n_right+n_down2+n_mod]=[source]*n_mod
    return

fRe=24.00*4
Nu=4.8608 #or some constant, laminar case
#hc=Nu*kt/D #will also be constant
n_tsteps=100000
dt=.1 # s
# consider adaptive time steps see Vijayan eq. (4.99)
beam_cycle=240 # s
beam_on=60 # s
alpha=kt/(rho_0*cp) #J/smk kgm2/ss2mK    kgm/s3K * 1/kg/m3 * 1/J/kgK
# alpha=0 # used to test possible oddities due to longitudinal conduction

def dbyds_plus(i,T,s):
    if(i==len(T)-1):
        ds_ave=((s[1]-s[0])+(s[len(T)-1]-s[len(T)-2]))/2
        derivative=(T[0]-T[len(T)-1])/ds_ave
    else:
        derivative=(T[i+1]-T[i])/(s[i+1]-s[i])
    return derivative

def dTds_plus(T,ds):
    derivative=(np.roll(T,-1)-T)/ds
    return derivative

def dbyds_minus(i,T,s):
    if(i==0):
        ds_ave=((s[1]-s[0])+(s[len(T)-1]-s[len(T)-2]))/2
        derivative=(T[0]-T[len(T)-1])/ds_ave
    else:
        derivative=(T[i]-T[i-1])/(s[i]-s[i-1])
    return derivative

def dTds_minus(T,ds):
    derivative=(T-np.roll(T,1))/(ds*1.)
    return derivative

def d2byds2(i,T,s):
    if(i==0):
        ds_ave=((s[1]-s[0])+(s[len(T)-1]-s[len(T)-2]))/2
    elif(i==len(T)-1):
        ds_ave=((s[1]-s[0])+(s[len(T)-1]-s[len(T)-2]))/2
    else:
        ds_ave=(s_array[i+1]-s_array[i-1])/2
    second_derivative=(dbyds_plus(i,T,s)-dbyds_minus(i,T,s))/ds_ave
    return second_derivative

def d2Tds2(T,ds):
    second_derivative=(dTds_plus(T,ds)-dTds_minus(T,ds))/ds
    return second_derivative


def dbyds_avg(i,T,s):
    return (dbyds_plus(i,T,s)+dbyds_minus(i,T,s))/2

def dbyds_directional(w,i,T,s):
    if(w>0):
        derivative=dbyds_minus(i,T,s)
    elif(w<0):
        derivative=dbyds_plus(i,T,s)
    else:
        derivative=0
    return derivative

def dTds_directional(w,T,ds):
    if(w>0):
        derivative=dTds_minus(T,ds)
    elif(w<=0):
        derivative=dTds_plus(T,ds)
    return derivative

def Gimme_Gamma():
    # Gamma is a purely geometrical factor Vijayan Eq. (4.25)
    return sum(ds_array/A_array)

Gamma=Gimme_Gamma()

def friction_factor(Re):
    flam=np.divide(64,Re,out=np.zeros_like(Re),where=Re>0.0000003) # assume cicular tube correlation if not in HEX
    flam_hex=np.divide(fRe,Re,out=np.zeros_like(Re),where=Re>0.0000003)
    flam[0:n_hex]=flam_hex[0:n_hex]
    ftr=1.2036*np.float_power(Re,-0.416,np.zeros_like(Re),where=Re>0.0000003) #from vijayan
    fturb=0.316*np.float_power(Re,-0.25,np.zeros_like(Re),where=Re>0.0000003)
    f=np.maximum(np.maximum(flam,ftr),fturb)
    f[Re<=0.0000003]=.01 # use a small amount of friction for smallest Re
    return f

def Nusselt(Re):
    Pr=mu*cp/kt
    # Gnielinski correlation https://en.wikipedia.org/wiki/Nusselt_number
    fgniel=(0.79*np.log(Re,np.zeros_like(Re),where=Re>0.00000003)-1.64)**-2
    Nugniel=(fgniel/8)*(Re-1000)*Pr/(1+12.7*(fgniel/8)**.5*(Pr**(2/3)-1))
    number=np.zeros_like(Re)
    number[:n_hex]=4.8608*(Re<3000)[:n_hex]
    number[n_hex:]=3.657*(Re<3000)[n_hex:]
    number=number+Nugniel*(Re>=3000)
    return number

for tstep in range(0,n_tsteps):
    t=dt*tstep
    # set beam current
    if(t>beam_cycle*1):
        if(t%beam_cycle<beam_on):
            set_beam_current(10.)
        else:
            set_beam_current(10.)
    else:
        set_beam_current(0.)

    # update fluid temperatures
    dTemp=dt*(-(w/rho_0)/A_array*dTds_minus(T_array,ds_array)+source_array+source_wall_to_fluid+alpha*d2Tds2(T_array,ds_array))
    T_array=T_array+dTemp

    # update wall temperatures
    dTemp_wall=dt*(source_fluid_to_wall+alpha_wall*d2Tds2(T_wall,ds_array))
    T_wall=T_wall+dTemp_wall


    # perform the integral of rho*dz
    rho_integral=np.sum(-rho_0*beta_t*g*T_array*(z_array-np.roll(z_array,1)))

    # calculate sum of f*L/(D_h*A) around the loop
    D_h_array=4*A_array/P_array
    Re_array=abs(4*w/(P_array*mu))
    f_array=friction_factor(Re_array)
    foa2=f_array*ds_array/(D_h_array*A_array**2)
    foa2_sum=np.sum(foa2)
    foa2K_sum=foa2_sum+KoA2cont1+KoA2cont2+2*KoA245+KoA2exp1+KoA2cont1+KoA290+KoA2valve+KoA2exp2

    friction_pressure_drops=foa2*w**2/(2*rho_0) # Pa?
    friction_term=foa2K_sum*w**2/(2*rho_0)

    # dw for this time step
    dw=(dt/Gamma)*(-friction_term-rho_integral) # Vijayan (4.25)
    w=w+dw

    # update heat transfer correlations

    hcw=Nusselt(Re_array)*kt/D_h_array # Barron eq'n 6.15

    source_wall_to_fluid=-4*hcw*(T_array-T_wall)/(D_h_array*rho_0*cp)
    source_wall_to_fluid[0:n_hex]=source_wall_to_fluid[0:n_hex]*(1-perimeter/P)
    
    source_fluid_to_wall=4*hcw*(T_array-T_wall)/(D_h_array*rhocp_wall)
    source_fluid_to_wall[0:n_hex]=source_fluid_to_wall[0:n_hex]*(1-perimeter/P)

    # update hex cold source
    source_array[0:n_hex]=-4*hcw[0:n_hex]*(T_array[0:n_hex]-T_cold)/(D_hex*rho_0*cp)*perimeter/P

    sparse=1000 # sparseness of standard output
    if(tstep%sparse==0):
        print('This is time %f and w is %f'%(t,w))
        print(min(T_array),max(T_array),foa2_sum)
        print([np.sum(friction_pressure_drops[i:i+n_per]) for i in range(0,len(friction_pressure_drops),n_per)])
        Tminvalue.append(min(T_array))
        Tmaxvalue.append(max(T_array))
        wvalue.append(w)
        tvalue.append(t)


savetxt('w.txt', wvalue, delimiter=',')

print('foa2_sum = %f m^-2'%foa2_sum)
print('Ts = %f K'%T_cold)
print('P = %f m^-1'%(4*hc*A/(D*w*cp)))
print('top of hex height = %f m'%top_z_hex)
print('bot of hex height = %f m'%(top_z_hex-L_hex))
print('top of mod = %f m'%(very_bottom_z+L_mod))
print('very bottom of mod = %f m'%very_bottom_z)
h=(top_z_hex+top_z_hex-L_hex)/2-(very_bottom_z+L_mod+very_bottom_z)/2
print('height approx = %f m'%h)
print('w est. = %f kg/s'%(2*rho_0**2*beta_t*g*q_mod_total*h/(foa2_sum*cp))**(1./3.))

plt.plot(tvalue,wvalue,'r:')
plt.ylabel('w (kg/s)')
plt.title('Mass Flux as a Function of Time t')
plt.xlabel('Time (s)')
plt.show()

plt.plot(tvalue,Tminvalue,'b:',label='Min Temp')
plt.plot(tvalue,Tmaxvalue,'r:',label='Max Temp')
plt.legend(loc='upper left')
plt.title('Temperature as a Function of Time')
plt.ylabel('Temperature (K)')
plt.xlabel('Time (s)')
plt.show()



plt.plot(T_array,'b:')
#plt.plot(T1,'r:')
plt.title('Temperature as a Function of Position of Index')
plt.ylabel('Temperature (K)')
plt.xlabel('Position')
plt.show()

plt.plot(T_wall,'b:')
#plt.plot(T1,'r:')
plt.title('Temperature of Wall as a Function of Position of Index')
plt.ylabel('Temperature (K)')
plt.xlabel('Position')
plt.show()


plt.plot(s_array,T_array,'b:')
#plt.plot(T1,'r:')
plt.title('Temperature as a Function of Physical Position Around the Loop')
plt.ylabel('Temperature (K)')
plt.xlabel('Position (m)')
#plt.show()

plt.plot(z_array,T_array,'b:')
#plt.plot(T1,'r:')
plt.title('Temperature as a Function of Height Around the Loop')
plt.ylabel('Temperature (K)')
plt.xlabel('Height Position (m)')
#plt.show()















