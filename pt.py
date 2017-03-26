import numpy as np
import os

# given redshift and flag return the 3rd PT results
def get_3rd_data(z, flag):        #flag = '8600', '8610', '8411'
	command = '/Users/dalongcheng/Dropbox/Projects/Code/Codes/3rd-PT/V-2.0/3rd_PT %s %.6e >./Result/%s_3rd_%.3f.txt'%(flag, z, flag, z)
	os.system(command)
	file_name = './Result/%s_3rd_%.3f.txt'%(flag, z)
	theory = np.loadtxt(file_name)
	command = 'rm ./Result/%s_3rd_%.3f.txt'%(flag, z)
	os.system(command)            #the pt output is deleted
	return theory
