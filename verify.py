import numpy as np
import matplotlib.pyplot as plt

import cosmo
import scaling
import pt

# It is designed to replace the highest redshift data with PT and scale it according to the data fitted alpha.
# here we concern only 8610 and 8600, which are most important.

def show_legend(ax, location='best'):
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=location)

def scale_wf(i, data, theory, file_alpha, flag_theory):
	if flag_theory == 0:
		wf_ref = scaling.wf(data[-1])                           # the reference wf of the highest redshift of data
	else:
		wf_ref = scaling.wf(theory)                             # use theoretical value, if flag = 1
	if i == len(data)-1:
		return wf_ref                                           # do no scaling for the last one
	else:
		ds     = file_alpha[i,0]                                # the gr at as
		df     = file_alpha[i,1]
		alpha  = file_alpha[i,2]
		wf_res = np.power(wf_ref, np.power(df/ds, alpha))
		return wf_res

def verify_wf(data, rz, flag):                                          # flag == '8600' and '8610' only, '8411' will be verified at other place
	file_in = './Result/%s_das_daf_alpha.txt'%(flag)
	file_alpha = np.loadtxt(file_in)
	theory = pt.get_3rd_data(rz[-1], flag)
	#do plot	
	fig = plt.figure(figsize=(8, 12))
	#the first is a-alpha relation	
	ax1 = fig.add_subplot(3, 2, 1)
	ax1.set_xlim(0.1, 1.)
	ax1.set_ylim(-2.02, -1.66)
	ax1.set_xlabel('a')
	ax1.set_ylabel(r'$-\alpha$')
	ax1.plot(file_alpha[:,1], -file_alpha[:,2], "bo-", label=flag)
	#the second one tests high redshift theory with data	
	ax2 = fig.add_subplot(3, 2, 2)
	ax2.set_xscale('log')
	ax2.set_xlabel('k[h/Mpc]')	
	ax2.set_ylabel(r"$\tilde{W}(k)$")
	ax2.set_ylim(0., 1.1)
	for i in np.arange(-1, -5, -1):
		temp_theory = pt.get_3rd_data(rz[i], flag)
		ax2.plot(data[i][:,0], scaling.wf(data[i]), "g-")
		ax2.plot(temp_theory[:,0], scaling.wf(temp_theory), "r--")
	#the scaling tests according to the data
	ax3 = fig.add_subplot(3, 2, 3)		
	ax3.set_xscale('log')
	ax3.set_xlabel('k[h/Mpc]')	
	ax3.set_ylabel(r"$\tilde{W}(k)$")
	ax3.set_ylim(0., 1.1)
	for i in range(len(data)):
		ax3.plot(data[i][:,0], scaling.wf(data[i]), "g-")
		ax3.plot(data[i][:,0], scale_wf(i, data, theory, file_alpha, 0), "r--")
	#the relative difference
	ax4 = fig.add_subplot(3, 2, 4)
	ax4.set_xscale('log')
	ax4.set_xlabel("k[h/Mpc]")
	ax4.set_ylabel("Predicted/Measured")
	ax4.set_ylim(0.95, 1.05)
	ax4.set_yticks(np.arange(0.95, 1.05, 0.02))	
	for i in range(len(data)):
		ax4.plot(data[i][:,0], scale_wf(i, data, theory, file_alpha, 0)/scaling.wf(data[i]), "g-")
	#the scaling tests according to the pt 
	ax5 = fig.add_subplot(3, 2, 5)		
	ax5.set_xscale('log')
	ax5.set_xlabel('k[h/Mpc]')	
	ax5.set_ylabel(r"$\tilde{W}(k)$")
	ax5.set_ylim(0., 1.1)
	for i in range(len(data)):
		ax5.plot(data[i][:,0], scaling.wf(data[i]), "g-")
		ax5.plot(data[i][:,0], scale_wf(i, data, theory, file_alpha, 1), "r--")	
	#the relative difference 
	ax6 = fig.add_subplot(3, 2, 6)
	ax6.set_xscale('log')
	ax6.set_xlabel("k[h/Mpc]")
	ax6.set_ylabel("Predicted/Measured")
	ax6.set_ylim(0.95, 1.05)
	ax6.set_yticks(np.arange(0.95, 1.05, 0.02))		
	for i in range(len(data)):
		ax6.plot(data[i][:,0], scale_wf(i, data, theory, file_alpha, 1)/scaling.wf(data[i]), "g-")
	show_legend(ax1)		
	fig.set_tight_layout(True)
	return fig                                              # return the fig




