import numpy as np
from scipy import optimize

import vals
import cosmo
import pt

#this gives the P_dt/P_dd, which is the window function
def wf(data):
	data_eff = data[:,2]/data[:,1]
	return data_eff

#this gives the eta
def eta(data):
	data_eff = data[:,3]*data[:,1]/np.power(data[:,2], 2) - 1.0
	return data_eff
	
def get_alpha(i, kmin, kmax, data, a):                     # kmax is the maximum k used.
	res = optimize.minimize(chi_sq, 2, (i, kmin, kmax, data, a), method='nelder-mead', options={'disp':False, 'xtol':1e-5, 'ftol':1e-8})
	return res.x

def chi_sq(alpha, i, kmin, kmax, data, da):                # i=0 for z=0, da is the growth factor
	cond1 = data[i][:,0] <= kmax
	cond2 = data[i][:,0] >= kmin
	cond = np.logical_and(cond1, cond2)
	data_1 = data[i][cond,:]
	data_2 = data[-1][cond,:]                              # the largest redshift one of data
	wf1 = wf(data_1)
	wf2 = wf(data_2)
	wf1_eff = np.power(wf1, np.power(da[i],  -alpha))
	wf2_eff = np.power(wf2, np.power(da[-1], -alpha))      # the largest redshfit one of data
	chi2 = np.power(wf1_eff-wf2_eff, 2.).sum()
	return chi2

def get_a_alpha(data, rz, obj_cosmo, flag):   
	da = obj_cosmo.gr(rz)                                  # the growth factor	
	das_daf_alpha = np.zeros((rz.shape[0]-1, 3))           # das is the growth factor of the highest redshift
	for i in range(das_daf_alpha.shape[0]):
		das_daf_alpha[i,0] = da[-1]
		das_daf_alpha[i,1] = da[i]
		das_daf_alpha[i,2] = get_alpha(i, 0., 0.5, data, da)
	file_out = './Result/%s_das_daf_alpha.txt'%(flag)
	np.savetxt(file_out, das_daf_alpha, fmt='%.6e', delimiter='\t')
