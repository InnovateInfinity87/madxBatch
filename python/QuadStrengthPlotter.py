# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:17:03 2016

@author: wvandepo
"""
import matplotlib.pyplot as plt
import numpy as np
import Constants as c


kqf1_start =      1.46516e-02
kqd_start =     -1.46338e-02
kqf1_end =       1.46572e-02
kqd_end =     -1.46349e-02

kqf1 = kqf1_start 
kqd  = kqd_start 

m_f = (kqf1_start - kqf1_end) / (1 - 34095)
m_d = (kqd_start -  kqd_end ) / (1 - 34095)

n_f = kqf1_start - m_f
n_d = kqd_start - m_d

plt.figure()

x=[]
y=[]
ypert=[]
ypert2=[]
ypert3=[]
for i in range(34095):
    x.append(i)
    y.append(m_f * i + n_f)
    if(i>5000):
        pert=-1e-5*np.sin(2*np.pi*i/500)*np.exp((5000-i)/500)
        pert2=-1e-5*np.sin(2*np.pi*i/500)*np.exp((5000-i)/1000)
        pert3=-2e-5*np.sin(2*np.pi*i/500)*np.exp((5000-i)/500)
    else:
        pert=0
        pert2=0
        pert3=0
    ypert.append(m_f * i + n_f + pert)
    ypert2.append(m_f * i + n_f + pert2)
    ypert3.append(m_f * i + n_f + pert3)



plt.plot(x,ypert3,label=r"Damped Sine: $\frac{\Delta I}{I}=5.6e-4$, 500 turns period, 500 turns damping")
plt.plot(x,ypert2,label=r"Damped Sine: $\frac{\Delta I}{I}=2.8e-4$, 500 turns period, 1000 turns damping")
plt.plot(x,ypert,label=r"Damped Sine: $\frac{\Delta I}{I}=2.8e-4$, 500 turns period, 500 turns damping")

plt.plot(x,y,label="Lineair tune sweep")
plt.legend(loc=4)
plt.ylabel("Quadrupole strength")
plt.xlabel("Turns")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()
#
#data = np.loadtxt("../twiss/tune_var.prt",skiprows=8,unpack=True)
#
#f, axarr = plt.subplots(2, 1,sharex=True, sharey=False)
#f.suptitle(r"Quadrupole strength in function of turns", fontsize=14)
#
#
#axarr[0].plot(data[0],data[1], label='KQF')
#axarr[0].set_ylabel(r'KQF')
#axarr[0].set_xlabel(r'Turns')
#
#axarr[1].plot(data[0],data[3], label='KQD')
#axarr[1].set_ylabel(r'KQD')
#axarr[1].set_xlabel(r'Turns')

#plt.legend()



#plt.xlim([0,500])

