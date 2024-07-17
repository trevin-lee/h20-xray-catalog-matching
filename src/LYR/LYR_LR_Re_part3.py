import numpy as np
import matplotlib.pyplot as plt

import LYR_functions as lyr
from LYR_input_parameters import *
from LYR_distrib_part2 import *



############################################################# LR and Reliability ###
print('\n\nLR and Reliability:')

neigh_idx_lr = tree.query_ball_point(r_out_tree, r=r_lr)

d=0
for i in range(len(neigh_idx_lr)):
	for j in range(len(neigh_idx_lr[i])):
		rain_tmp = r_in_tree[neigh_idx_lr[i],0]
		decin_tmp = r_in_tree[neigh_idx_lr[i],1]
		Dra_tmp = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_tmp = (r_out_tree[i,1] - decin_tmp[j])
		rr_lr = lyr.quadratic_sum(Dra_tmp, Ddec_tmp)
		if rr_lr <= r_lr:
			d=d+1
	prog = 100*i/len(ID_output)
	print ('Counting sources: '+str(int(prog)+1)+'%', end="\r")

lr_inRA = np.zeros(d)
lr_inDEC = np.zeros(d)
lr_outRA = np.zeros(d)
lr_outDEC = np.zeros(d)
lr_mag_in = np.zeros(d)
lr_LR = np.zeros(d)
lr_inID = np.zeros(d)
lr_outID = np.zeros(d)
lr_r = np.zeros(d)
lr_flux_output = np.zeros(d)
lr_fluxerr_output = np.zeros(d)
lr_flux_input = np.zeros(d)
lr_fluxerr_input = np.zeros(d)
rel = []
lr_noID = []
flag = []
fr_r_tmp = np.zeros(d)
fr_fr_tmp = np.zeros(d)
fdf_delta_tmp = np.zeros(d)
fdf_fdf_tmp = np.ones(d)
nnm_ = np.zeros(d)
qm_ = np.zeros(d)


#fig99, ax99 = plt.subplots()
#rel_sum = []
d=0
for i in range(len(neigh_idx_lr)):
	lr_singleXID = []
	for j in range(len(neigh_idx_lr[i])):
		rain_tmp = r_in_tree[neigh_idx_lr[i],0]
		decin_tmp = r_in_tree[neigh_idx_lr[i],1]
		Dra_tmp = (r_out_tree[i,0] - rain_tmp[j])
		Ddec_tmp = (r_out_tree[i,1] - decin_tmp[j])
		rr_lr = lyr.quadratic_sum(Dra_tmp, Ddec_tmp)
		if rr_lr <= r_lr:
			ra_tmp = ra_input[neigh_idx_lr[i]]
			lr_inRA[d] = ra_tmp[j]
			dec_tmp = dec_input[neigh_idx_lr[i]]
			lr_inDEC[d] = dec_tmp[j]
			mag_tmp = mag_input[neigh_idx_lr[i]]
			lr_mag_in[d] = mag_tmp[j]
			inID_tmp = ID_input[neigh_idx_lr[i]]
			lr_inID[d] = inID_tmp[j]
			lr_outID[d] = ID_output[i]
			lr_outRA[d] = ra_output[i]
			lr_outDEC[d] = dec_output[i]
			lr_r[d] = rr_lr


			#calcolo la LR(fr,nm,qm) per ogni sorgente
			if delta_f2 == True:

				lr_flux_output[d] = np.power(10,flux_output[i])
				lr_fluxerr_output[d] = np.power(10,fluxerr_output[i])
				lr_flux_input[d] = np.power(10,lr_mag_in[d])

				if delta_f2_fluxerr == 0:
					if magerr_input == 0:
						lr_fluxerr_input[d] = 0
					elif magerr_input != 0:
						lr_fluxerr_input[d] = np.power(10,magerr_input)
				elif delta_f2_fluerr == 1:
					lr_fluxerr_input[d] = np.power(10,magerr_input)

				if (input_cat_type == 2) or (input_cat_type == 0):
					fr_r_tmp[d], fr_fr_tmp[d] = lyr.FR(sigma_output[i],sigma_input,Dra_tmp,Ddec_tmp)
					fdf_delta_tmp[d], fdf_fdf_tmp[d] = lyr.Delta_flux(lr_flux_output[d],lr_flux_input[d],lr_fluxerr_output[d],lr_fluxerr_input[d])
					nnm_[d] = nnm_interp(lr_mag_in[d])
					qm_[d] = qm_interp(lr_mag_in[d])
					lr_LR[d] = lyr.LR_(fr_fr_tmp[d],nnm_[d],qm_[d],fdf_fdf_tmp[d])
					#lr_LR[d] = lr_LR[d]*pnos[d]

				elif input_cat_type == 1:
					magerr_tmp = magerr_input[neigh_idx_lr[i]]
					lr_fluxerr_input[d] = np.power(10,magerr_tmp[j])
					sigma_in_tmp = sigma_input[neigh_idx_lr[i]]
					fr_r_tmp[d], fr_fr_tmp[d] = lyr.FR(sigma_output[i],sigma_in_tmp[j],Dra_tmp,Ddec_tmp)
					fdf_delta_tmp[d], fdf_fdf_tmp[d] = lyr.Delta_flux(lr_flux_output[d],lr_flux_input[d],lr_fluxerr_output[d],lr_fluxerr_input[d])
					nnm_[d] = nnm_interp(lr_mag_in[d])
					qm_[d] = qm_interp(lr_mag_in[d])
					lr_LR[d] = lyr.LR_(fr_fr_tmp[d],nnm_[d],qm_[d],fdf_fdf_tmp[d])

				# only for the outuput file
				lr_flux_output[d] = np.log10(lr_flux_output[d])
				lr_fluxerr_output[d] = np.log10(lr_fluxerr_output[d])


			elif delta_f2 == False:

				lr_flux_output[d] = 0.
				lr_fluxerr_output[d] = 0.
				lr_flux_input[d] = 0.

				if (input_cat_type == 2) or (input_cat_type == 0):
					fr_r_tmp[d], fr_fr_tmp[d] = lyr.FR(sigma_output[i],sigma_input,Dra_tmp,Ddec_tmp)
					nnm_[d] = nnm_interp(lr_mag_in[d])
					qm_[d] = qm_interp(lr_mag_in[d])
					lr_LR[d] = lyr.LR(fr_fr_tmp[d],nnm_[d],qm_[d])
				elif input_cat_type == 1:
					sigma_in_tmp = sigma_input[neigh_idx_lr[i]]
					fr_r_tmp[d], fr_fr_tmp[d] = lyr.FR(sigma_output[i],sigma_in_tmp[j],Dra_tmp,Ddec_tmp)
					nnm_[d] = nnm_interp(lr_mag_in[d])
					qm_[d] = qm_interp(lr_mag_in[d])
					lr_LR[d] = lyr.LR(fr_fr_tmp[d],nnm_[d],qm_[d])

			else:
				print('\nProblem with delta_f2 option!\n')

			#... e memorizzo quelle della singola XID temporaneamente:
			lr_singleXID = np.append(lr_singleXID, lr_LR[d])
			d=d+1

	#calcolo Re con le LR nell'intorno delle singole XID:
	if len(lr_singleXID) != 0.:
		summLR = np.sum(lr_singleXID)
		single_outID_LR = []
		rel_sum_tmp = []
		for k in range(len(lr_singleXID)):
			rel_tmp = lyr.Re(summLR, lr_singleXID[k], Q)
			rel = np.append(rel, rel_tmp)
			#rel_sum_tmp = np.append(rel_sum_tmp, rel_tmp)
			single_outID_LR = np.append(single_outID_LR, lr_singleXID[k])
		# find the maximum likelihood of each output source (1 == max LR; 0 == non max LR):
		LR_max_tmp = max(single_outID_LR)
		#rel_sum = np.append(rel_sum, np.sum(rel_sum_tmp))
		for ff in range(len(single_outID_LR)):
			if single_outID_LR[ff] < LR_max_tmp:
				flag_tmp = 0
				flag = np.append(flag, flag_tmp)
			elif single_outID_LR[ff] == LR_max_tmp:
				flag_tmp = 1
				flag = np.append(flag, flag_tmp)
	else:
		lr_noID = np.append(lr_noID, ID_output[i])
	prog = 100*i/len(ID_output)
	print ('Filling array: '+str(int(prog)+1)+'%', end="\r")
print()


fig5, ax = plt.subplots()
ax.scatter(lr_inRA, lr_inDEC, color='b', marker='o', alpha=0.5, label='Mathced input sources')
ax.scatter(ra_output, dec_output, color='r', marker='.', label='All Output sources')
#ax.scatter(lr_outRA, lr_outDEC, color='gold', marker='.', s=2, label='Output sources matched')
ax.set_xlabel('RA [deg]')
ax.set_ylabel('DEC [deg]')
ax.legend()
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
plt.gca().invert_xaxis()

if save_images == True:
	fig5.savefig(path_images+'match_results_'+str(int(r_in_tm))+'rlr'+str(int(r_lr))+add_str+'.png',
				 bbox_inches="tight", dpi=250)

non_matched=len(lr_noID)/len(ID_output)
print('Unmatched',len(lr_noID),'sources (',int(100*non_matched)+1,'%):\n', lr_noID)


if save_output == True:
	outinfo = open(path_output+filename_outinfo,"a")
	outinfo.write('\n' + str(len(lr_noID)) + ' unmatched' + ' sources ('
				  + str(int(100*non_matched)+1) + ' %): \n' + str(lr_noID) + '\n')
	outinfo.close()

	width=12
	out_file = open(path_output+filename_LR,"w")
	out_file.write("#" + '{0: <{width}.8s}'.format("OutID", width=width)
					+ '{0: <{width}.8s}'.format("RAout", width=width)
					+ '{0: <{width}.8s}'.format("DECout", width=width)
					+ '{0: <{width}.8s}'.format("mag_out", width=width)
					+ '{0: <{width}.8s}'.format("InID", width=width)
					+ '{0: <{width}.8s}'.format("RAin", width=width)
					+ '{0: <{width}.8s}'.format("DECin", width=width)
					+ '{0: <{width}.8s}'.format("mag_in", width=width)
					+ '{0: <{width}.8s}'.format("r", width=width)
					+ '{0: <{width}.8s}'.format("f(r)", width=width)
					+ '{0: <{width}.8s}'.format("df", width=width)
					+ '{0: <{width}.8s}'.format("f(df)", width=width)
					+ '{0: <{width}.8s}'.format("q(m)", width=width)
					+ '{0: <{width}.8s}'.format("n(m)", width=width)
					+ '{0: <{width}.8s}'.format("LR", width=width)
					+ '{0: <{width}.8s}'.format("Rel", width=width)
					+ '{0: <{width}.8s}'.format("flag", width=width)
					+ "\n")
	out_file.close()

	out_file = open(path_output+filename_LR,"a")
	for d in range(len(lr_inID)):
		out_file.write('{0: <{width}.0f}'.format(lr_outID[d], width=width+1)
					   + '{0: <{width}.6f}'.format(lr_outRA[d], width=width)
					   + '{0: <{width}.6f}'.format(lr_outDEC[d], width=width)
					   + '{0: <{width}.4f}'.format(lr_flux_output[d], width=width)
					   + '{0: <{width}.0f}'.format(lr_inID[d], width=width)
					   + '{0: <{width}.6f}'.format(lr_inRA[d], width=width)
					   + '{0: <{width}.6f}'.format(lr_inDEC[d], width=width)
					   + '{0: <{width}.4f}'.format(lr_mag_in[d], width=width)
					   + '{0: <{width}.3f}'.format(fr_r_tmp[d], width=width)
					   + '{0: <{width}.3e}'.format(fr_fr_tmp[d], width=width)
					   + '{0: <{width}.3e}'.format(fdf_delta_tmp[d], width=width)
					   + '{0: <{width}.3e}'.format(fdf_fdf_tmp[d], width=width)
					   + '{0: <{width}.3e}'.format(qm_[d], width=width)
					   + '{0: <{width}.3e}'.format(nnm_[d], width=width)
					   + '{0: <{width}.3e}'.format(lr_LR[d], width=width)
					   + '{0: <{width}.3f}'.format(rel[d], width=width)
					   + '{0: <{width}.0f}'.format(flag[d], width=width)
					   + "\n")
	out_file.close()
	print('... done.\n')


#plt.show()
