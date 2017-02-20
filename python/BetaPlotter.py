# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 16:20:54 2016

@author: wvandepo
"""

import matplotlib.pyplot as plt
#import numpy as np
from MADXreader import TwissRead

f, axarr = plt.subplots(4, 1,sharex=True, sharey=False)
f.suptitle(r"Optical characteristics before and after thinning", fontsize=14)

variable_data1, parameters1 = TwissRead('twiss_before_thinning.prt')
variable_data2, parameters2 = TwissRead('twiss_after_thinning.prt')

axarr[0].plot(variable_data1['S'],variable_data1['BETX'], label='Before')
axarr[0].plot(variable_data2['S'],variable_data2['BETX'], label='After')
axarr[0].set_ylabel(r'$\beta_x$')
axarr[0].set_xlabel(r'$S$[m]')
axarr[0].axhline(0, color='black')

axarr[1].plot(variable_data1['S'],variable_data1['BETY'], label='Before')
axarr[1].plot(variable_data2['S'],variable_data2['BETY'], label='After')
axarr[1].set_ylabel(r'$\beta_y$')
axarr[1].set_xlabel(r'$S$[m]')
axarr[1].axhline(0, color='black')

axarr[2].plot(variable_data1['S'],variable_data1['X'], label='Before')
axarr[2].plot(variable_data2['S'],variable_data2['X'], label='After')
axarr[2].set_ylabel(r'$x$')
axarr[2].set_xlabel(r'$S$[m]')
axarr[2].axhline(0, color='black')

axarr[3].plot(variable_data1['S'],variable_data1['Y'], label='Before')
axarr[3].plot(variable_data2['S'],variable_data2['Y'], label='After')
axarr[3].set_ylabel(r'$y$')
axarr[3].set_xlabel(r'$S$[m]')
axarr[3].axhline(0, color='black')

plt.xlim([1500,2000])
plt.legend()
