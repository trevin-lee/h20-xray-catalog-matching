import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.spatial import KDTree

import LYR_functions as lyr
from LYR_input_parameters import *
from LYR_sigma_part1 import *

# outuput file info:
if save_output == True:
	outinfo = open(path_output+filename_outinfo,"a")
	outinfo.write('\nSearching radii:\n'
				  + 'r_fr = ' + str(r_fr) + ';  r_in_tn = ' + str(r_in_tm)
				  + ';  r_min_nm = ' +str(r_min_nm) + ';  r_max_nm = ' +str(r_max_nm)
				  + ';  r_lr = ' +str(r_lr) + '\n')
	outinfo.close()

###################################################################### f(r) ########
fig1, ax = plt.subplots()
print('\nFinding FR...')

r_in_tree = np.column_stack((ra_input*3600, dec_input*3600))
r_out_tree = np.column_stack((ra_output*3600, dec_output*3600))

tree = KDTree(r_in_tree)

neighbors_idx = tree.query_ball_point(r_out_tree, r=r_fr)

for i in range(len(neighbors_idx)):
	for j in range(len(neighbors_idx[i])):
		rain_tmp = r_in_tree[neighbors_idx[i],0]
		decin_tmp = r_in_tree[neighbors_idx[i],1]
		Dra_fr = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_fr = (r_out_tree[i,1] - decin_tmp[j])
		rr_fr = lyr.quadratic_sum(Dra_fr, Ddec_fr)
		sigma_out_fr = sigma_output[i]
		if (input_cat_type == 2) or (input_cat_type == 0):
			sigma_in_fr = sigma_input
		elif (input_cat_type == 1):
			sigma_in_fr = sigma_input[j]

		# call FR:
		r_tot, fr_tot = lyr.FR(sigma_out_fr, sigma_in_fr, Dra_fr, Ddec_fr)
		ax.scatter(r_tot, fr_tot, c='r', marker='.', zorder=3)

	prog = 100*i/len(neighbors_idx)
	print (str(int(prog)+1)+'%', end="\r")

ax.set_yscale('log')
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.set_xlabel('r')
ax.set_ylabel('f(r)')
ax.grid(ls=':', color='grey', alpha=0.4, zorder=0)

print('\n... done.')


###################################################################### n(m) ########
print('\nSearching n(m)...')

neigh_idx_nmax = tree.query_ball_point(r_out_tree, r=r_max_nm)

a=0
for i in range(len(neigh_idx_nmax)):
	for j in range(len(neigh_idx_nmax[i])):
		rain_tmp = r_in_tree[neigh_idx_nmax[i],0]
		decin_tmp = r_in_tree[neigh_idx_nmax[i],1]
		Dra_nm = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_nm = (r_out_tree[i,1] - decin_tmp[j])
		rr_nm = lyr.quadratic_sum(Dra_nm, Ddec_nm)
		if rr_nm > r_min_nm:
			a=a+1
	prog = 100*i/len(neigh_idx_nmax)
	print ('Counting sources: ' + str(int(prog)+1) + '%', end="\r")

nn_mag_in = np.zeros(a)
nn_in_ID = np.zeros(a)
a=0
for i in range(len(neigh_idx_nmax)):
	for j in range(len(neigh_idx_nmax[i])):
		rain_tmp = r_in_tree[neigh_idx_nmax[i],0]
		decin_tmp = r_in_tree[neigh_idx_nmax[i],1]
		Dra_nm = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_nm = (r_out_tree[i,1] - decin_tmp[j])
		rr_nm = lyr.quadratic_sum(Dra_nm, Ddec_nm)
		if rr_nm > r_min_nm:
			tmp_ID = ID_input[neigh_idx_nmax[i]]
			nn_in_ID[a] = tmp_ID[j]
			tmp_mag = mag_input[neigh_idx_nmax[i]]
			nn_mag_in[a] = tmp_mag[j]
			a=a+1

	prog = 100*i/len(neigh_idx_nmax)
	print ('Filling arrays ' + str(int(prog)+1) + '%', end="\r")


unique_n, indices = np.unique(nn_in_ID, return_index=True)
n_mag_in = np.zeros(len(indices))
unique_n2 = np.zeros(len(indices))
for i in range(len(indices)):
	n_mag_in[i] = nn_mag_in[indices[i]]
	unique_n2[i] = nn_in_ID[indices[i]]

print('\n... done.')


###################################################################### q(m) ########
print('\nSearching q(m)...')
acircle = lyr.acircle(r_in_tm)
aannu = lyr.aarea(r_min_nm, r_max_nm)

#-------------------------------------------------------------- total(m):

neigh_idx_totm = tree.query_ball_point(r_out_tree, r=r_in_tm)

a=0
for i in range(len(neigh_idx_totm)):
	for j in range(len(neigh_idx_totm[i])):
		rain_tmp = r_in_tree[neigh_idx_totm[i],0]
		decin_tmp = r_in_tree[neigh_idx_totm[i],1]
		Dra_tm = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_tm = (r_out_tree[i,1] - decin_tmp[j])
		rr_tm = lyr.quadratic_sum(Dra_tm, Ddec_tm)
		if rr_tm <= r_in_tm:
			a=a+1

	prog = 100*i/len(neigh_idx_totm)
	print ('Counting sources: ' + str(int(prog)+1) + '%', end="\r")

tot_mag_in = np.zeros(a)
tot_in_ID = np.zeros(a)
a=c=d=0
for i in range(len(neigh_idx_totm)):
	b=0
	for j in range(len(neigh_idx_totm[i])):
		rain_tmp = r_in_tree[neigh_idx_totm[i],0]
		decin_tmp = r_in_tree[neigh_idx_totm[i],1]
		Dra_tm = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_tm = (r_out_tree[i,1] - decin_tmp[j])
		rr_tm = lyr.quadratic_sum(Dra_tm, Ddec_tm)
		if rr_tm <= r_in_tm:
			ID_tmp = ID_input[neigh_idx_totm[i]]
			tot_in_ID[a] = ID_tmp[j]
			magin_tmp = mag_input[neigh_idx_totm[i]]
			tot_mag_in[a] = magin_tmp[j]
			a=a+1
			b=b+1
	if b != 0:
		c=c+1
	elif b == 0:
		d=d+1
	prog = 100*i/len(neigh_idx_totm)
	print ('Filling arrays: ' + str(int(prog)+1) + '%', end="\r")


unique_tot, indices = np.unique(tot_in_ID, return_index=True)
total_mag_in = np.zeros(len(indices))
unique_tot2 = np.zeros(len(indices))
for i in range(len(indices)):
	total_mag_in[i] = tot_mag_in[indices[i]]
	unique_tot2[i] = tot_in_ID[indices[i]]

print('\n... done.')


r0 = r_in_tm
N1 = c
Nx = len(ID_output)
#Q = N1/Nx
Q = 1 - (d/Nx)
print('\nNin,tot=',N1, '; Nout=', Nx, '; Q=', round(Q,2))

if save_output == True:
	outinfo = open(path_output+filename_outinfo,"a")
	outinfo.write('\n' + 'Nin,tot = ' + str(N1) + ';  Nout = ' + str(Nx) + ';  Q = ' +str(round(Q,2)) + '\n')
	outinfo.close()



fig3, ax = plt.subplots()
bb=0
for b in range(len(mag_input)):
	if mag_input[b] != nomag:
		bb=bb+1
tmp_mag_input = np.zeros(bb)
bb=0
for b in range(len(mag_input)):
	if mag_input[b] != nomag:
		tmp_mag_input[bb] = mag_input[b]
		bb=bb+1
bins = np.linspace(min(tmp_mag_input),max(tmp_mag_input),distrib_bins)
bar_width = bins[1]-bins[0]
bin_height, bin_edges, patches = ax.hist(n_mag_in, bins=bins, color='r', alpha=0.5)
ax.set_ylabel('num. sources')
ax.set_xlabel('input magnitude/flux')
ax.set_title('Total background n(m)')
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)


# ... total(m):
fig4, ax = plt.subplots()
area_ratio = aannu/acircle
#rescaling n(m) area to tot(m) area:
bin_height = np.array(bin_height/area_ratio, dtype=float)
bins_ = bins + (bins[1]-bins[0])/2
bins_ = np.delete(bins_,-1)

bin_h, bin_e, patches = ax.hist(total_mag_in, bins=bins, density=0, histtype='step',
								linestyle='--', linewidth=1.5, color='dodgerblue', label='total(m)')
x_ = bins_
y_ = bin_h
#smoothing:
Xplot_totm = np.linspace(min(tmp_mag_input),max(tmp_mag_input),1000)[:, np.newaxis]
total_mag_in = total_mag_in[:, np.newaxis]
# la kde funziona su distribuzioni normalizzate:
kde = KernelDensity(kernel='gaussian', bandwidth=bar_width).fit(total_mag_in)
log_dens = kde.score_samples(Xplot_totm)
Yplot_totm = np.exp(log_dens)*sum(bin_h)*bar_width

# background n(m):
y_nm,x_nm,patches = ax.hist(lyr.bartohist(bins_, bin_height), color='mediumblue', density=0, histtype='step',
							linestyle='-.', linewidth=1.5, bins=bins, label='bkg n(m)')

x_ = bins_
y_ = bin_height
Xplot_nm = Xplot_totm
n_mag_in = n_mag_in[:, np.newaxis]
kde = KernelDensity(kernel='gaussian', bandwidth=bar_width).fit(n_mag_in)
log_dens = kde.score_samples(Xplot_nm)
Yplot_nm = np.exp(log_dens)*sum(bin_height)*bar_width
#ax.plot(Xplot_nm[:, 0], Yplot_nm, ':', color='mediumblue', label='n(m)')
nm_interp = interpolate.interp1d(Xplot_nm[:, 0], Yplot_nm)


# real(m):
Xplot_realm = Xplot_totm
Yplot_realm = Yplot_totm - Yplot_nm
real_height = bin_h - bin_height
real_height = np.array(np.around(real_height,decimals=0), dtype=int)
for i in range(len(real_height)):
	if real_height[i] < 0.:
		real_height[i] = 0.

real_hhist = np.around(bin_h,0) - np.around(bin_height,0)
real_xhist = lyr.bartohist(bins_,real_hhist)
y_realm,x_realm,patches = ax.hist(real_xhist, color='orange', density=0, histtype='stepfilled', alpha=0.5, bins=bins, label='real(m)')
real_xhist = real_xhist[:, np.newaxis]
Xplot_realm = np.linspace(min(tmp_mag_input),max(tmp_mag_input),1000)[:, np.newaxis]
kde = KernelDensity(kernel='gaussian', bandwidth=bar_width).fit(real_xhist)
log_dens = kde.score_samples(Xplot_realm)
Yplot_realm = np.exp(log_dens)*sum(y_realm)*bar_width
ax.plot(Xplot_realm[:, 0], Yplot_realm, '-', color='k', label='real(m)')


ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
ax.set_ylabel('num. sources')
ax.set_xlabel('Input magnitude/flux')
ax.set_title('Counterparts magnitude distribution' + add_title)
ax.legend(loc='upper left')
#ax.grid(ls=':', color='grey', alpha=0.5)
#plt.gca().invert_xaxis()
if save_images == True:
	fig4.savefig(path_images+'cparts_mag_distrib_rt'+str(int(r_in_tm))+'rlr'+str(int(r_lr))+add_str+'.png',
				 bbox_inches="tight", dpi=250)




# q(m) = real(m)/somm(real(m))*Q normalizzo real(m) al numero di sorgenti
summ_real = np.sum(y_realm)
qm = (Yplot_realm/summ_real)*round(Q,2)
#qm = Yplot_realm*round(Q,1)
Xplot_qm = Xplot_realm[:, 0]
qm_interp = interpolate.interp1d(Xplot_qm, qm)

# stessa cosa per n(m):
nnm = Yplot_nm/(np.pi*np.power(r0,2)*Nx)
Xplot_nnm = Xplot_nm[:, 0]
nnm_interp = interpolate.interp1d(Xplot_nnm, nnm)


#plt.show()
