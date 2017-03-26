import cosmo
import vals

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# Given a sigma^2(M,z) using the last redshift 5.09 of 8600 as reference to calculate a strating neff. 
# Then find the corresponding z to neff of 8610. Then use the transformation function to obtain a new
# index-redshift relation for 8610. We do a universal test of the results by plotting three figures:
# (1) index over growth factor,
# (2) index over neff of the final redshift, the starting neff is different
# (3) index over neff of the finial redshift,
# (4) index under the limitation of same starting and final redshift.

# for interpolate.InterpolatedUnivariateSpline, the input x array must be increasing.

def init_states():
	data_8600 = np.loadtxt('./Result/8600_das_daf_alpha.txt')
	data_8610 = np.loadtxt('./Result/8610_das_daf_alpha.txt')	
	vals.Sp_alpha_a_8600 = interpolate.InterpolatedUnivariateSpline(data_8600[::-1,1], data_8600[::-1,2])
	vals.Sp_alpha_a_8610 = interpolate.InterpolatedUnivariateSpline(data_8610[::-1,1], data_8610[::-1,2])

def beta(x, z2, rz, Sp):
	a1 = 1.0/(1.0+rz[-1])
	a2 = 1.0/(1.0+z2)
	y = Sp(x) * (np.log10(x) - np.log10(a1))
	y += Sp(a2)*np.log10(a1/a2)
	y /= np.log10(x/a2)
	return y

def beta_limit(x, z2, rz, Sp):	
	a1 = 1.0/(1.0+rz[-1])	
	a2 = 1.0/(1.0+z2)
	y = Sp(x) + x*np.log(x/a1)*Sp(x, nu=1)
	return y	

# this function calculate \beta(a2, a2)
def universal_test(sigma2, obj_cosmo_8600, obj_cosmo_8610):
	init_states()
	obj_cosmo_8600.setup_spline_da_neff(sigma2, '8600')
	obj_cosmo_8610.setup_spline_da_neff(sigma2, '8610')
	astart_8600 = 1.0 /(1.0+vals.z_8600[-1])
	zstart_8600 = vals.z_8600[-1]
	zneff_8600  = obj_cosmo_8600.Sp_neff_da(astart_8600)   # the starting zneff of 8600
	anew_8610   = obj_cosmo_8610.Sp_da_neff(zneff_8600)    # the new starting a of 8610
	znew_8610   = 1.0/anew_8610 - 1.0
	da_8600 = obj_cosmo_8600.gr(vals.z_8600[:-1])          # omit the last one
	da_8610 = obj_cosmo_8610.gr(vals.z_8610[:-1])
	astart_8610 = 1.0/(1.0+vals.z_8610[-1])
	zneff_8610 = obj_cosmo_8610.Sp_neff_da(astart_8610)
	#do some plots
	fig = plt.figure(figsize=(10, 10))	
	ax1 = fig.add_subplot(2, 2, 1)
	ax1.set_xlabel("a")
	ax1.set_ylabel(r"$\alpha$")
	x1 = da_8600
	y1 = vals.Sp_alpha_a_8600(x1)
	x2 = da_8610
	y2 = vals.Sp_alpha_a_8610(x2)
	ax1.plot(da_8600, y1, "ro-")
	ax1.plot(da_8610, y2, "bo-")
	ax2 = fig.add_subplot(2, 2, 2)
	ax2.set_xlabel(r"$n_{eff}$")
	ax2.set_ylabel(r"$\alpha$")
	ax2.set_title(r"$\sigma^2=%.2f, n^{I}_{8600}=%.2f, n^{I}_{8610}=%.2f$"%(sigma2, zneff_8600, zneff_8610))
	x1 = obj_cosmo_8600.Sp_neff_da(da_8600)
	y1 = vals.Sp_alpha_a_8600(da_8600)
	x2 = obj_cosmo_8610.Sp_neff_da(da_8610)
	y2 = vals.Sp_alpha_a_8610(da_8610)
	ax2.plot(x1, y1, "ro-")
	ax2.plot(x2, y2, "bo-")
	ax3 = fig.add_subplot(2, 2, 3)
	ax3.set_xlabel(r"$n_{eff}$")
	ax3.set_ylabel(r"$\alpha$")
	ax3.set_title(r"$\sigma^2=%.2f, n^{I}_{eff}=%.2f, Z^{I}_{8610}=%.2f$"%(sigma2, zneff_8600, znew_8610))
	x1 = obj_cosmo_8600.Sp_neff_da(da_8600)
	y1 = vals.Sp_alpha_a_8600(da_8600)
	x2 = obj_cosmo_8610.Sp_neff_da(da_8610)
	y2 = beta(da_8610, znew_8610, vals.z_8610, vals.Sp_alpha_a_8610)
	ax3.plot(x1, y1, "ro-")
	ax3.plot(x2, y2, "bo-")
	ax4 = fig.add_subplot(2, 2, 4)
	ax4.set_xlabel(r"$n_{eff}$")
	ax4.set_ylabel(r"$\lim\beta$")
	x1 = obj_cosmo_8600.Sp_neff_da(da_8600)
	y1 = beta_limit(da_8600, zstart_8600, vals.z_8600, vals.Sp_alpha_a_8600)
	x2 = obj_cosmo_8610.Sp_neff_da(da_8610)
	y2 = beta_limit(da_8610, znew_8610, vals.z_8610, vals.Sp_alpha_a_8610)
	ax4.plot(x1, y1, "ro-")
	ax4.plot(x2, y2, "bo-")
	fig.set_tight_layout(True)
	return fig
