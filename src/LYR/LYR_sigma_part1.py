import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import LYR_functions as lyr
from LYR_input_parameters import *


##################################### output file info:
if save_output == True:
	outinfo = open(path_output+filename_outinfo,"w")
	outinfo.write('###################\n'+'### OUTPUT INFO ###\n' + '###################\n\n')
	outinfo.close()

#################################################################### std deviations ########
print('\nCalculating std...')
if delta_f2 == True:
	fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(10,4))
else:
	fig, (ax1,ax2) = plt.subplots(1,2, figsize=(10,4))


############################################################################# output(X-ray):
err_output = np.sqrt(np.power(ra_err_output_s,2)+np.power(dec_err_output_s,2))
err_output = err_output*plus_erf*3600

entries, bin_edges, patches = ax1.hist(err_output, bins=20, color='dodgerblue', histtype='stepfilled',
									   alpha=0.7, label='Output errors', zorder=2)
bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
# gaussian fit method 1
parameters1, cov_matrix1 = curve_fit(lyr.gaussian, bin_middles, entries)
xgauss = np.linspace(min(err_output), max(err_output), 1000)
ygauss = lyr.gaussian(xgauss, *parameters1)
ax1.plot(xgauss, ygauss, ls='--', lw=1, color='blue')
label='$\sigma_{out}$='+str(abs(round(parameters1[1],2)))

ax1.set_title('Output: ' + label)
ax1.grid(c='grey', ls=':', alpha=0.6, zorder=0)
ax1.set_xlabel('[arcsec]')
ax1.legend(fontsize=9)

msigma_output = round(parameters1[1],2)

if sigma_out_finder == '1sigma':
	sigma_output = np.zeros(len(err_output))
	for i in range(len(err_output)):
		sigma_output[i] = msigma_output
	print('Output mean sigma:', abs(round(msigma_output,2)))
elif sigma_out_finder == 'all':
	sigma_output = err_output
	print('Using all output sigma.')



############################################################################### input(ONIR):
if input_cat_type == 2:
	x_rahist = []
	x_dechist = []

	for i in range(len(ID1_input_s)):
		if (ra2_input_s[i] != noxy) and (dec2_input_s[i] != noxy) and (mag2_input_s[i] != nomag):
			if (ra1_input_s[i] != noxy) and (dec1_input_s[i] != noxy) and (mag1_input_s[i] != nomag):
				hist_Dra = (ra1_input_s[i]*3600 - ra2_input_s[i]*3600)
				hist_Ddec = (dec1_input_s[i]*3600 - dec2_input_s[i]*3600)
				x_rahist = np.append(x_rahist, hist_Dra)
				x_dechist = np.append(x_dechist, hist_Ddec)
	label1 = '$\Delta$RA'
	label2 = '$\Delta$DEC'

elif input_cat_type == 1:
	x_hist=[]
	for i in range(len(ID1_input_s)):
			if (ra1_input_s[i] != noxy) and (dec1_input_s[i] != noxy) and (mag1_input_s[i] != nomag):
				err_input = lyr.quadratic_sum(ra1_input_s, dec1_input_s)
				x_hist = np.append(err_input*3600, hist_tmp)
	label='Input errors'

elif input_cat_type == 0:
	x_hist = 0
	label='Input errors'



if input_cat_type == 2:
	bins = np.linspace(min(min(x_rahist),min(x_dechist)), max(max(x_rahist),max(x_dechist)), 20)

	ax2.hist(x_rahist, bins=bins, color='darkorange', alpha=0.3, zorder=3)
	entries, bin_edges, patches = ax2.hist(x_rahist, bins=bins, color='darkorange', histtype='step', label=label1, zorder=3)
	bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
	parameters, cov_matrix = curve_fit(lyr.gaussian, bin_middles, entries)
	xgauss = np.linspace(min(x_rahist), max(x_rahist), 1000)
	ygauss = lyr.gaussian(xgauss, *parameters)
	title1 = ' $\sigma_{RA,in}$='+str(abs(round(parameters[1],2)))
	ax2.plot(xgauss, ygauss, ls='--', color='r', zorder=4)
	sigma_rain = abs(parameters[1])

	ax2.hist(x_dechist, bins=bins, color='orchid', alpha=0.3, zorder=3)
	entries, bin_edges, patches = ax2.hist(x_dechist, bins=bins, color='orchid', histtype='step', label=label2, zorder=3)
	bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
	parameters, cov_matrix = curve_fit(lyr.gaussian, bin_middles, entries)
	xgauss = np.linspace(min(x_dechist), max(x_dechist), 1000)
	ygauss = lyr.gaussian(xgauss, *parameters)
	title2 = ' $\sigma_{DEC,in}$='+str(abs(round(parameters[1],2)))
	ax2.plot(xgauss, ygauss, ls='--', color='purple', lw=1, zorder=4)
	sigma_decin = abs(parameters[1])

	sigma_input = lyr.quadratic_sum(sigma_rain, sigma_decin)
	title3 = '$\sigma_{r,in}$=' + str(abs(round(sigma_input,2)))
	print('Input mean sigma:', round(sigma_input,2))

elif input_cat_type == 1:
	sigma_input = lyr.quadratic_sum(ra1_input_s_err,dec1_input_s_err)
	print('Using all input sigma errors.')

elif input_cat_type == 0:
	sigma_input = sigma_input
	title0 = '$\sigma_{r,in}$=' + str(abs(round(sigma_input,2)))
	print('Input sigma:', round(sigma_input,2))



if input_cat_type == 2:
	ax2.set_title('Input: ' + title1 + ' ' + title2 + ' ' + title3)
elif input_cat_type == 1:
	ax2.set_title('Input')
elif input_cat_type == 0:
	ax2.set_title('Input ' + title0)
ax2.grid(c='grey', ls=':', alpha=0.6, zorder=0)
ax2.legend(fontsize=9)
ax2.set_xlabel('[arcsec]')




######################################################################### flux output distrib:
if delta_f2 == True:
	bins = np.linspace(min(flux_output_s), max(flux_output_s), 20)
	entries, bin_edges, patches = ax3.hist(flux_output_s, bins=bins, color='forestgreen', histtype='stepfilled',
										   alpha=0.7, label='Output flux errors', zorder=2)
	bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
	# gaussian fit method 1
	parameters1, cov_matrix1 = curve_fit(lyr.gaussian, bin_middles, entries)
	xgauss = np.linspace(min(flux_output_s), max(flux_output_s), 1000)
	ygauss = lyr.gaussian(xgauss, *parameters1)
	ax3.plot(xgauss, ygauss, ls='--', color='g')
	ax3.plot([], [], color='white', label='$\sigma_{f,out}$='+str(abs(round(parameters1[1],2))))

	ax3.set_title("Output")
	ax3.grid(c='grey', ls=':', alpha=0.6, zorder=0)
	ax3.set_xlabel('log f')
	ax3.legend(fontsize=9, loc='upper right')

	mfsigma_output = round(parameters1[1],2)
	print('Output flux mean sigma:', abs(round(mfsigma_output,2)))



if save_images == True:
	fig.savefig(path_images+'sigma_inout_'+str(int(r_in_tm))+'rlr'+str(int(r_lr))+add_str+'.png',
				bbox_inches="tight", dpi=250)

print('... done.\n')
