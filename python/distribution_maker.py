# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:07:27 2016

@author: wvandepo
"""
import make_ps_distribution as dis

dis.get_gauss_distribution(output='../input/initial_distribution_gauss.txt', input='twiss_after_thinning.prt',n_part=1000, sigmas=4, beam_t='FT', seed=254)