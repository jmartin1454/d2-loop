#!/usr/bin/python

from math import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy import interpolate
from scipy.misc import derivative
from numpy import asarray
from numpy import savetxt
from scipy.optimize import curve_fit


#Values of Q_h for Vij 4.20

q=np.array([10,20,30,40,50,60,70,80,90,100])

#values of w found by running Vij 4.20
ws=np.array([log10(0.003016475538632849),log10(0.0038005210276165022),log10(0.004350510549436035),log10(0.0047883564432621255),log10(0.005158100614800805),log10(0.005481299819024171),log10(0.005770310119920915),log10(0.006032951077265698),log10(0.006274521970541968),log10(0.006498799542063221)])

#Values of Q_h for Vij 4.21
qT=(2/3)*np.array([log10(10),log10(20),log10(30),log10(40),log10(50),log10(60),log10(70),log10(80),log10(90),log10(100)])


#values of T found by running Vij 4.21
T=np.array([0.5049698651598858,0.8015896951670426,1.0503796476477798,1.2724443253550317,1.4765408429904652,1.6673737576420722,1.8478391109691001,2.019879460639543,2.1848777131351054,2.3438624874370766])



#Values of w found by running my program for diff q_mod_total
dws=np.array([0.003013,0.003946,0.004612,0.005149,0.005607,0.006009,0.006371,0.006701,0.007006,0.007290])
#Values of T found by running my program for diff q_mod_total
dTs=(-1)*np.array([19.810336525302503-20.314501232180135,19.829673876645472-20.600431327640923,19.85483517963164-20.84425158366479,19.884088984306917-21.066026274158595,19.91645590261472-21.273605062367853,19.951312454499302-21.471038506144833,19.988229370602983-21.660742447309485,20.026894524000067-21.844296338836866,20.067071588860973-22.022800600693394,20.10857596884922-22.19705871351428])


#plt.plot(ws,q,'r:')
#plt.yscale('log')
#plt.xscale('log')
#plt.ylabel('w')
##plt.title('log(w) vs 1/3*log(q)')
#plt.xlabel('q')
#plt.grid(True)
#plt.show()
##
#plt.plot(T,q,'r:')
#plt.yscale('log')
#plt.xscale('log')
#plt.ylabel('dT')
#plt.grid(True)
##plt.title('log(dT) vs 2/3*log(q)')
#plt.xlabel('q')
#plt.show()

def powerlaw(x,c,m):
    return c*x**m

plt.plot(q,dws,'ro')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('$w$ (kg/s)')
plt.grid(True)
#plt.
#plt.ylim(1e-4,1e-1)
#plt.xlim(1,1000)
#plt.title('log(w) vs 1/3*log(q)')
plt.xlabel('$q$ (W)')
popt,pcov=curve_fit(powerlaw,q,dws)
print(popt)
plt.plot(q,powerlaw(q,*popt),'--',label='$w=%6.4f q^{%3.2f}$'%tuple(popt))
plt.legend()
plt.show()
#
plt.plot(q,dTs,'ro')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('$\Delta T$ (K)')
plt.grid(True)
#plt.title('log(dT) vs 2/3*log(q)')
plt.xlabel('$q$ (W)')
popt,pcov=curve_fit(powerlaw,q,dTs)
print(popt)
plt.plot(q,powerlaw(q,*popt),'--',label='$\Delta T=%6.4f q^{%3.2f}$'%tuple(popt))
plt.legend()
plt.show()

