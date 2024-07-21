from matplotlib import pyplot as plt
import numpy as np
import time

from sklearn.neighbors import KernelDensity
from scipy import interpolate
from scipy.spatial import KDTree
from scipy.optimize import curve_fit

import sys; sys.path.append("../lyr")
import LYR_functions as lyr

class LYR:

    def __init__():
        pass

    def run_lyr(
            r_in_tm,
            r_fr,
            r_min_nm, 
            r_max_nm,
            data_input,
            data_output,
            data1,
            distrib_bins,
        ):

        r_lr = r_fr
        plus_erf = 1
        input_cat_type = 0
        nomag = -99.
        ID_input = np.array(data_input[:,0])
        ra_input = np.array(data_input[:,1])
        dec_input = np.array(data_input[:,2])
        mag_input = np.array(data_input[:,3])
        ID_output = np.array(data_output[:,0])
        ra_output = np.array(data_output[:,1])
        dec_output = np.array(data_output[:,2])
        ra_err_output = np.array(data_output[:,3]) / 3600
        dec_err_output = np.array(data_output[:,4]) / 3600
        sigma_input = 0 # 0.3 per cosmos
        ra_err_output_s = ra_err_output
        dec_err_output_s = dec_err_output
        toc = time.time()
        err_output = np.sqrt(np.power(ra_err_output_s,2)+np.power(dec_err_output_s,2))
        err_output = err_output*plus_erf*3600
        entries, bin_edges = np.histogram(err_output, bins=20)
        bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
        parameters1, cov_matrix1 = curve_fit(lyr.gaussian, bin_middles, entries)
        xgauss = np.linspace(min(err_output), max(err_output), 1000)
        ygauss = lyr.gaussian(xgauss, *parameters1)
        msigma_output = round(parameters1[1],2)
        sigma_output = np.zeros(len(err_output))
        for i in range(len(err_output)):
            sigma_output[i] = msigma_output
        r_in_tree = np.column_stack((ra_input*3600, dec_input*3600))
        r_out_tree = np.column_stack((ra_output*3600, dec_output*3600))
        tree = KDTree(r_in_tree)
        neighbors_idx = tree.query_ball_point(r_out_tree, r=r_fr, workers=8)
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

            prog = 100*i/len(neighbors_idx)
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
        unique_n, indices = np.unique(nn_in_ID, return_index=True)
        n_mag_in = np.zeros(len(indices))
        unique_n2 = np.zeros(len(indices))
        for i in range(len(indices)):
            n_mag_in[i] = nn_mag_in[indices[i]]
            unique_n2[i] = nn_in_ID[indices[i]]



        ###################################################################### q(m) ########
        acircle = lyr.acircle(r_in_tm)
        aannu = lyr.aarea(r_min_nm, r_max_nm)

        #-------------------------------------------------------------- total(m):

        neigh_idx_totm = tree.query_ball_point(r_out_tree, r=r_in_tm, workers=8)

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


        unique_tot, indices = np.unique(tot_in_ID, return_index=True)
        total_mag_in = np.zeros(len(indices))
        unique_tot2 = np.zeros(len(indices))
        for i in range(len(indices)):
            total_mag_in[i] = tot_mag_in[indices[i]]
            unique_tot2[i] = tot_in_ID[indices[i]]


        r0 = r_in_tm
        N1 = c
        Nx = len(ID_output)
        #Q = N1/Nx
        Q = 1 - (d/Nx)

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
        bin_height, bin_edges = np.histogram(n_mag_in, bins=bins)

        area_ratio = aannu/acircle
        #rescaling n(m) area to tot(m) area:
        bin_height = np.array(bin_height/area_ratio, dtype=float)
        bins_ = bins + (bins[1]-bins[0])/2
        bins_ = np.delete(bins_,-1)

        bin_h, bin_e = np.histogram(total_mag_in, bins=bins, density=0)
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
        y_nm, x_nm = np.histogram(lyr.bartohist(bins_, bin_height), density=0, bins=bins)

        x_ = bins_
        y_ = bin_height
        Xplot_nm = Xplot_totm
        n_mag_in = n_mag_in[:, np.newaxis]
        kde = KernelDensity(kernel='gaussian', bandwidth=bar_width).fit(n_mag_in)
        log_dens = kde.score_samples(Xplot_nm)
        Yplot_nm = np.exp(log_dens)*sum(bin_height)*bar_width
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
        y_realm, x_realm = np.histogram(real_xhist, density=0, bins=bins)
        real_xhist = real_xhist[:, np.newaxis]
        Xplot_realm = np.linspace(min(tmp_mag_input),max(tmp_mag_input),1000)[:, np.newaxis]
        kde = KernelDensity(kernel='gaussian', bandwidth=bar_width).fit(real_xhist)
        log_dens = kde.score_samples(Xplot_realm)
        Yplot_realm = np.exp(log_dens)*sum(y_realm)*bar_width

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


        ############################################################# LR and Reliability ###
        neigh_idx_lr = tree.query_ball_point(r_out_tree, r=r_lr, workers=8)

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
        rel = []
        lr_noID = []
        flag = []
        fr_r_tmp = np.zeros(d)
        fr_fr_tmp = np.zeros(d)
        nnm_ = np.zeros(d)
        qm_ = np.zeros(d)

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

                    #... e memorizzo quelle della singola XID temporaneamente:
                    lr_singleXID = np.append(lr_singleXID, lr_LR[d])
                    d=d+1

            #calcolo Re con le LR nell'intorno delle singole XID:
            if len(lr_singleXID) != 0.:
                summLR = np.sum(lr_singleXID)
                single_outID_LR = []
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


        ############################################################# Completness ###

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


        a=b=0
        for i in range(len(lr_LR)):
            if lr_LR[i] > LRth[RCmax_index]:
                a=a+1
                if flag[i] == 1:
                    b=b+1

        tic = time.time()
        tt = tic - toc
        tt_str = 'seconds.'
        if tt > 60:
            tt = tt/60
            tt_str = 'minutes.'
        elif tt > 3600:
            tt = tt/3600
            tt_str = 'hours.'

        LRth = LRth[RCmax_index]
        CLRth = C[RCmax_index]
        RLRth = R[RCmax_index]

        return LRth, CLRth, RLRth
