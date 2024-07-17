import matplotlib.pyplot as plt
import numpy as np
from LYR_input_parameters import *


if cursors == True:
	from mpldatacursor import datacursor

filename = path_output+filename_LR
do,di,rain,decin,magin,r,fr,df,fdf,qm,nm,lr,rel = np.genfromtxt(filename, usecols=(0,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)

fig, ax = plt.subplots(1,2, figsize=(10,4))

for i in range(len(lr)):
	ax[0].scatter(r[i], lr[i], marker='.', color='b')
	ax[1].scatter(r[i], rel[i], marker='.', color='r')

ax[0].set_yscale('log')
ax[0].set_xlabel('r [arcsec]')
ax[0].set_ylabel('LR')
ax[0].grid(ls=':', color='grey', alpha=0.4, zorder=0)
ax[1].set_xlabel('r [arcsec]')
ax[1].set_ylabel('Reliability')
ax[1].grid(ls=':', color='grey' ,alpha=0.4, zorder=0)
fig.suptitle('LR and reliability' + add_title)
ax[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)



if save_images == True:
	fig.savefig(path_images+'LRrel_'+str(int(r_in_tm))+'rlr'+str(int(r_lr))+add_str+'.png',
				bbox_inches="tight", dpi=250)

fig, ax = plt.subplots()
for i in range(len(lr)):

	di_ = int(di[i])
	rain_ = round(rain[i],6)
	decin_ = round(decin[i],6)
	do_ = int(do[i])
	magin_ = round(magin[i],2)
	r_ = round(r[i],2)
	fr_ = round(fr[i],3)
	df_ = round(df[i],2)
	fdf_ = round(fdf[i],3)
	nm_ = round(nm[i],3)
	qm_ = round(qm[i],3)
	lr_ = round(lr[i],3)
	rel_ = round(rel[i],3)

	if lr[i] != 0.:
		if cursors == True:
			ax.scatter(lr[i], rel[i], marker='.', color='mediumblue',
					  label='IDout: {}'.format(do_)+'\nIDint: {}'.format(di_)+'\nLR: {}'.format(lr_)+'\nRe: {}'.format(rel_)
					  + '\nq(m): {}'.format(qm_) + '\nn(m): {}'.format(nm_) + '\nf(r): {}'.format(fr_) + '\nf(df): {}'.format(fdf_)
					  + '\nr: {}'.format(r_) + '\ndf: {}'.format(df_) + '\nmag: {}'.format(magin_), zorder=3)
		elif cursors== False:
			ax.scatter(lr[i], rel[i], marker='.', color='mediumblue')



ax.set_xscale('log')
ax.set_xlabel('LR')
ax.set_ylabel('Reliability')
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.grid(ls=':', color='grey', alpha=0.4)
ax.set_title(add_title)

if cursors == True:
	datacursor(xytext=(15, 15), display='multiple', draggable=True,
			   bbox=dict(fc='white', alpha=1),formatter='{label}'.format)

if save_images == True:
	fig.savefig(path_images+'LRandrel_'+str(int(r_in_tm))+'rlr'+str(int(r_lr))+add_str+'.png',
				bbox_inches="tight", dpi=250)

#plt.show()
