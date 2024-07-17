import numpy as np
import matplotlib.pyplot as plt

import LYR_functions as lyr
from LYR_input_parameters import *
from LYR_LR_Re_part3 import *


############################################################# Completness ###
print('\nCompleteness:')

#LRth = np.arange(1e-10,1.,0.000001)
LRth = np.arange(0.01,100.,0.05)

RC = []
R = []
C = []

for i in range(len(LRth)):
	compl_LR = []
	compl_rel = []
	l=0
	for j in range(len(lr_LR)):
		if lr_LR[j] > LRth[i]:
			compl_LR = np.append(compl_LR, lr_LR[j])
			compl_rel = np.append(compl_rel, rel[j])
			l=l+1

	R_ = np.sum(compl_rel)/l
	R = np.append(R, R_)
	C_ = np.sum(compl_rel)/Nx
	C = np.append(C, C_)
	RC_ = (R_+C_)/2
	RC = np.append(RC, RC_)

RCmax_index = np.where(RC==np.nanmax(RC))

RCmax_index = RCmax_index[0]
if len(RCmax_index) > 1:
	RCmax_index = int(RCmax_index[0])
else:
	RCmax_index = int(RCmax_index)
print( '\nLRth: ' + str(round(LRth[RCmax_index],2)))
print( 'C[LRth]: ' + str(round(C[RCmax_index],2)))
print( 'R[LRth]: ' + str(round(R[RCmax_index],2)))

a=b=0
for i in range(len(lr_LR)):
	if lr_LR[i] > LRth[RCmax_index]:
		a=a+1
		if flag[i] == 1:
			b=b+1
print(str(a) + '/' + str(len(lr_LR)) + ' total input sources with LR > LRth (' + str(round((100*a/len(lr_LR)),2)) + '%)')
print(str(b) + '/' + str(Nx) + ' output sources with a counterpart with LR > LRth (' + str(round((100*b/Nx),2)) + '%)\n')

if save_output == True:
	outinfo = open(path_output+filename_outinfo,"a")
	outinfo.write('\n(R+C)/2 max: ' + str(round(np.max(RC),2)) + '\n'
				 + 'LRth: ' + str(round(LRth[RCmax_index],2)) + '\n'
				 + 'C[LRth]: ' + str(round(C[RCmax_index],2)) + '\n'
				 + 'R[LRth]: ' + str(round(R[RCmax_index],2)) + '\n'
				 + str(a) + '/' + str(len(lr_LR)) + ' total sources with LR > LRth (' + str(round((100*a/len(lr_LR)),2)) + '%)\n'
				 + str(b) + '/' + str(Nx) + ' output sources with a counterpart with LR > LRth (' + str(round((100*b/Nx),2)) + '%)\n')
	outinfo.close()

print('... done.\n')

fig6, ax = plt.subplots()

ax.plot(LRth, RC, color='r', label='(R+C)/2', zorder=3)
ax.plot(LRth, R, color='mediumblue', ls='--', label='R', zorder=3)
ax.plot(LRth, C, color='mediumblue', ls=':', label='C', zorder=3)
ax.axvline(color='r', lw=1, ls='--', alpha=0.6, x=LRth[RCmax_index])

ax.set_xlabel('$LR_{th}$')
ax.grid(ls=':', color='grey', alpha=0.3, zorder=0)
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.set_title('Completness' + add_title)
ax.set_xscale('log')
ax.legend()

if save_images == True:
	fig6.savefig(path_images+'completeness_'+str(int(r_in_tm))+'rlr'+str(int(r_lr))+add_str+'.png',
				 bbox_inches="tight", dpi=250)

