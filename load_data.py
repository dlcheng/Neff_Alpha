import numpy as np

import vals

def step_to_redshift(step):                          # works for 8600, 8610, and 8411
	total_step = 1200.
	init_z = 72.0
	dp = init_z/total_step
	j_factor = 1.0 + init_z
	pa = 1.0 + dp * step
	a = pa / j_factor
	z = 1.0/a - 1
	return z	

def load_data_set(flag):
	if flag == '8600':
		nb_24 = np.loadtxt("../Data/Power/Bin-20/wf8600_0130_z_7.30_2048.txt")		
		nb_23 = np.loadtxt("../Data/Power/Bin-20/wf8600_0183_z_5.09_2048.txt")
		nb_22 = np.loadtxt("../Data/Power/Bin-20/wf8600_0255_z_3.48_2048.txt")
		nb_21 = np.loadtxt("../Data/Power/Bin-20/wf8600_0354_z_2.28_2048.txt")
		nb_20 = np.loadtxt("../Data/Power/Bin-20/wf8600_0439_z_1.67_2048.txt")
		nb_19 = np.loadtxt("../Data/Power/Bin-20/wf8600_0463_z_1.54_2048.txt")
		nb_18 = np.loadtxt("../Data/Power/Bin-20/wf8600_0489_z_1.41_2048.txt")
		nb_17 = np.loadtxt("../Data/Power/Bin-20/wf8600_0516_z_1.28_2048.txt")
		nb_16 = np.loadtxt("../Data/Power/Bin-20/wf8600_0544_z_1.17_2048.txt")
		nb_15 = np.loadtxt("../Data/Power/Bin-20/wf8600_0574_z_1.06_2048.txt")
		nb_14 = np.loadtxt("../Data/Power/Bin-20/wf8600_0605_z_0.96_2048.txt")
		nb_13 = np.loadtxt("../Data/Power/Bin-20/wf8600_0638_z_0.86_2048.txt")
		nb_12 = np.loadtxt("../Data/Power/Bin-20/wf8600_0673_z_0.76_2048.txt")
		nb_11 = np.loadtxt("../Data/Power/Bin-20/wf8600_0709_z_0.68_2048.txt")
		nb_10 = np.loadtxt("../Data/Power/Bin-20/wf8600_0748_z_0.59_2048.txt")
		nb_09 = np.loadtxt("../Data/Power/Bin-20/wf8600_0788_z_0.51_2048.txt")
		nb_08 = np.loadtxt("../Data/Power/Bin-20/wf8600_0831_z_0.44_2048.txt")
		nb_07 = np.loadtxt("../Data/Power/Bin-20/wf8600_0876_z_0.36_2048.txt")
		nb_06 = np.loadtxt("../Data/Power/Bin-20/wf8600_0923_z_0.29_2048.txt")
		nb_05 = np.loadtxt("../Data/Power/Bin-20/wf8600_0973_z_0.23_2048.txt")
		nb_04 = np.loadtxt("../Data/Power/Bin-20/wf8600_1025_z_0.17_2048.txt")
		nb_03 = np.loadtxt("../Data/Power/Bin-20/wf8600_1080_z_0.11_2048.txt")
		nb_02 = np.loadtxt("../Data/Power/Bin-20/wf8600_1138_z_0.05_2048.txt")
		nb_01 = np.loadtxt("../Data/Power/Bin-20/wf8600_1200_z_0.00_2048.txt")
		step  = np.array([1200,1138,1080,1025,973,923,876,831,788,748,709,673,638,605,574,544,516,489,463,439,354,255,183])
		vals.data_8600  = [nb_01,nb_02,nb_03,nb_04,nb_05,nb_06,nb_07,nb_08,nb_09,nb_10,nb_11,nb_12,nb_13,nb_14,nb_15,nb_16,nb_17,nb_18,nb_19,nb_20,nb_21,nb_22,nb_23]
		vals.z_8600 = step_to_redshift(step)
	if flag == '8610':
		nb_24 = np.loadtxt("../Data/Power/Bin-20/wf8610_0130_z_7.30_2048.txt")		
		nb_23 = np.loadtxt("../Data/Power/Bin-20/wf8610_0183_z_5.09_2048.txt")		
		nb_22 = np.loadtxt("../Data/Power/Bin-20/wf8610_0255_z_3.48_2048.txt")
		nb_21 = np.loadtxt("../Data/Power/Bin-20/wf8610_0354_z_2.28_2048.txt")
		nb_20 = np.loadtxt("../Data/Power/Bin-20/wf8610_0439_z_1.67_2048.txt")
		nb_19 = np.loadtxt("../Data/Power/Bin-20/wf8610_0463_z_1.54_2048.txt")
		nb_18 = np.loadtxt("../Data/Power/Bin-20/wf8610_0489_z_1.41_2048.txt")
		nb_17 = np.loadtxt("../Data/Power/Bin-20/wf8610_0516_z_1.28_2048.txt")
		nb_16 = np.loadtxt("../Data/Power/Bin-20/wf8610_0544_z_1.17_2048.txt")
		nb_15 = np.loadtxt("../Data/Power/Bin-20/wf8610_0574_z_1.06_2048.txt")
		nb_14 = np.loadtxt("../Data/Power/Bin-20/wf8610_0605_z_0.96_2048.txt")
		nb_13 = np.loadtxt("../Data/Power/Bin-20/wf8610_0638_z_0.86_2048.txt")
		nb_12 = np.loadtxt("../Data/Power/Bin-20/wf8610_0673_z_0.76_2048.txt")
		nb_11 = np.loadtxt("../Data/Power/Bin-20/wf8610_0709_z_0.68_2048.txt")
		nb_10 = np.loadtxt("../Data/Power/Bin-20/wf8610_0748_z_0.59_2048.txt")
		nb_09 = np.loadtxt("../Data/Power/Bin-20/wf8610_0788_z_0.51_2048.txt")
		nb_08 = np.loadtxt("../Data/Power/Bin-20/wf8610_0831_z_0.44_2048.txt")
		nb_07 = np.loadtxt("../Data/Power/Bin-20/wf8610_0876_z_0.36_2048.txt")
		nb_06 = np.loadtxt("../Data/Power/Bin-20/wf8610_0923_z_0.29_2048.txt")
		nb_05 = np.loadtxt("../Data/Power/Bin-20/wf8610_0973_z_0.23_2048.txt")
		nb_04 = np.loadtxt("../Data/Power/Bin-20/wf8610_1025_z_0.17_2048.txt")
		nb_03 = np.loadtxt("../Data/Power/Bin-20/wf8610_1080_z_0.11_2048.txt")
		nb_02 = np.loadtxt("../Data/Power/Bin-20/wf8610_1138_z_0.05_2048.txt")
		nb_01 = np.loadtxt("../Data/Power/Bin-20/wf8610_1200_z_0.00_2048.txt")
		step  = np.array([1200,1138,1080,1025,973,923,876,831,788,748,709,673,638,605,574,544,516,489,463,439,354,255])
		vals.data_8610  = [nb_01,nb_02,nb_03,nb_04,nb_05,nb_06,nb_07,nb_08,nb_09,nb_10,nb_11,nb_12,nb_13,nb_14,nb_15,nb_16,nb_17,nb_18,nb_19,nb_20,nb_21,nb_22]
		vals.z_8610 = step_to_redshift(step)
	if flag == '8411':
		nb_01 = np.loadtxt("../Data/Power/Bin-20/wf8411_1200_z_0.00_2048.txt")
		nb_02 = np.loadtxt("../Data/Power/Bin-20/wf8411_0788_z_0.51_2048.txt")
		nb_03 = np.loadtxt("../Data/Power/Bin-20/wf8411_0605_z_0.96_2048.txt")
		step  = np.array([1200, 788, 605])
		vals.data_8411  = [nb_01, nb_02, nb_03]
		vals.z_8411 = step_to_redshift(step)

def load_all_data():
	load_data_set('8600')
	load_data_set('8610')
	load_data_set('8411')

