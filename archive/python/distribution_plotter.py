# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 13:46:53 2016

@author: wvandepo
"""
import make_ps_distribution as dis
import matplotlib.pyplot as plt
import itertools
from matplotlib import rcParams
rcParams.update({'font.size': 6})

#dis.get_fat_halo(output='initial_distribution', input="twiss_before_track_thin.tfs",n_halo=(4, 6), seed=2, beam_t='FT', n_part=10)
#dis.get_gauss_distribution(output='initial_distribution_Gauss', input='twiss_before_track_thin.tfs',sigmas=5, beam_t='FT', n_part=1000)


f, axarr = plt.subplots(4, 4,sharex=False, sharey=False)
f.suptitle("Initial distributions", fontsize=10)

##############################################################################
x, px, y, py, pt = dis.get_gauss_distribution(output='test', input='twiss_after_thinning.prt',n_part=1000, sigmas=4, beam_t='FT', seed=254)
i=0
axarr[i, 0].plot(x,y,".")
axarr[i, 0].set_title('Gauss sigma=4')
axarr[i, 1].plot(x,px,".")
axarr[i, 2].plot(y,py,".")

axarr[i,0].set_xlabel(r'$X$ [m]')
axarr[i,0].set_ylabel(r'$Y$ [m]')
axarr[i,1].set_xlabel(r'$X$ [m]')
axarr[i,1].set_ylabel(r'$p_X$ [GeV]')
axarr[i,2].set_xlabel(r'$Y$ [m]')
axarr[i,2].set_ylabel(r'$p_Y$ [GeV]')
axarr[i,0].grid(True)
axarr[i,1].grid(True)
axarr[i,2].grid(True)

axarr[i,3].hist(pt,50,normed=1,alpha=.75,histtype='stepfilled')
axarr[i,3].set_xlabel(r'$\Delta E$ [GeV]')
axarr[i,3].set_ylabel(r'$P(\Delta E)$')
axarr[i,3].grid(True)

###############################################################################
x, px, y, py, pt = dis.get_fat_halo(output='test', input='twiss_after_thinning.prt', n_part=1000, n_halo=(4, 6), beam_t='FT', seed=257861)
i=1
axarr[i, 0].plot(x,y,".")
axarr[i, 0].set_title('Fat halo, n_halo=(4, 6)')
axarr[i, 1].plot(x,px,".")
axarr[i, 2].plot(y,py,".")

axarr[i,0].set_xlabel(r'$X$')
axarr[i,0].set_ylabel(r'$Y$')
axarr[i,1].set_xlabel(r'$X$')
axarr[i,1].set_ylabel(r'$p_X$')
axarr[i,2].set_xlabel(r'$Y$')
axarr[i,2].set_ylabel(r'$p_Y$')
axarr[i,0].grid(True)
axarr[i,1].grid(True)
axarr[i,2].grid(True)

axarr[i,3].hist(pt,50,normed=1,alpha=.75,histtype='stepfilled')
axarr[i,3].set_xlabel(r'$\Delta E$ [GeV]')
axarr[i,3].set_ylabel(r'$P(\Delta E)$')
axarr[i,3].grid(True)

###############################################################################
x, px, y, py, pt = dis.get_halo_distribution(output='test', input='twiss_after_thinning.prt', n_part=1000, n_halo=4, beam_t='FT', seed=2742)
i=2
axarr[i, 0].plot(x,y,".")
axarr[i, 0].set_title('Halo, n_halo=4')
axarr[i, 1].plot(x,px,".")
axarr[i, 2].plot(y,py,".")

axarr[i,0].set_xlabel(r'$X$')
axarr[i,0].set_ylabel(r'$Y$')
axarr[i,1].set_xlabel(r'$X$')
axarr[i,1].set_ylabel(r'$p_X$')
axarr[i,2].set_xlabel(r'$Y$')
axarr[i,2].set_ylabel(r'$p_Y$')
axarr[i,0].grid(True)
axarr[i,1].grid(True)
axarr[i,2].grid(True)

axarr[i,3].hist(pt,50,normed=1,alpha=.75,histtype='stepfilled')
axarr[i,3].set_xlabel(r'$\Delta E$ [GeV]')
axarr[i,3].set_ylabel(r'$P(\Delta E)$')
axarr[i,3].grid(True)

###############################################################################
x, px, y, py, pt = dis.get_halo_distribution(output='test', input='twiss_after_thinning.prt', n_part=1000, n_halo=6, beam_t='FT', seed=225469)
i=3
axarr[i, 0].plot(x,y,".")
axarr[i, 0].set_title('Halo, n_halo=6')
axarr[i, 1].plot(x,px,".")
axarr[i, 2].plot(y,py,".")

axarr[i,0].set_xlabel(r'$X$')
axarr[i,0].set_ylabel(r'$Y$')
axarr[i,1].set_xlabel(r'$X$')
axarr[i,1].set_ylabel(r'$p_X$')
axarr[i,2].set_xlabel(r'$Y$')
axarr[i,2].set_ylabel(r'$p_Y$')
axarr[i,0].grid(True)
axarr[i,1].grid(True)
axarr[i,2].grid(True)

axarr[i,3].hist(pt,50,normed=1,alpha=.75,histtype='stepfilled')
axarr[i,3].set_xlabel(r'$\Delta E$ [GeV]')
axarr[i,3].set_ylabel(r'$P(\Delta E)$')
axarr[i,3].grid(True)

# Turn on the proper x or y axes ticks.
for i in range(4):
    for j in range(4):
        for tick in axarr[i,j].get_xticklabels():
            tick.set_rotation(45)
        
plt.tight_layout()
plt.show()