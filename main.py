# This routine can calculate the effective index N_eff defined as
# N_eff = 3(-dlog\sigma^2(M, z)/dlogM - 1)
# The definition is consistent for scale free case of either Gaussian and Top-hat 
# window function. Here we implement only the Top-hat window for general power spectrum.
#                                                        - dlcheng. 2015/10/20
import constants
import cosmo
import load_data
import vals
import scaling
import verify
import pt
import universal_test

import numpy as np
import matplotlib.pyplot as plt

reload(constants)
reload(cosmo)
reload(load_data)
reload(scaling)
reload(verify)
reload(pt)
reload(universal_test)


#set parameters for 8600, 8610, 8411
cos_8600 = cosmo.cosmo(1.0,   0.0,   0.166, 1.0,   0.83, 2, 1e2, 1e19, 1000)
cos_8610 = cosmo.cosmo(1.0,   0.0,   0.166, 0.968, 0.83, 4, 1e2, 1e19, 1000)
cos_8411 = cosmo.cosmo(0.268, 0.732, 0.166, 0.968, 0.83, 4, 1e2, 1e19, 1000)

#load data
load_data.load_all_data()

#do the scaling fit of data
scaling.get_a_alpha(vals.data_8600, vals.z_8600, cos_8600, '8600')
scaling.get_a_alpha(vals.data_8610, vals.z_8610, cos_8610, '8610')
scaling.get_a_alpha(vals.data_8411, vals.z_8411, cos_8411, '8411')

#verify the data fit results
fig_8600 = verify.verify_wf(vals.data_8600, vals.z_8600, '8600')
fig_8610 = verify.verify_wf(vals.data_8610, vals.z_8610, '8610')

#tests of universality
fig_universal = universal_test.universal_test(0.2, cos_8600, cos_8610)


