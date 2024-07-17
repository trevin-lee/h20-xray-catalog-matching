import numpy as np
import os
import LYR_functions as lyr
import matplotlib.pyplot as plt


#from Hickox+06 (teta : offaxis)
def PSF_radius(teta):
	r90 = 1 + 10*np.power(((teta*3600)/600),2)  # in arcsec
	return r90/3600  # in degrees

#from Puccetti+09
def Xerr(r_psf, cts):
	err = r_psf/np.sqrt(cts)
	return err


#################################### off-axis calc:
# on-axis (J1030+0524) in degrees
ra_center = 157.6130002
dec_center = 5.4153190


file = "/Users/Alessandro/Documents/post_TESI/codes/LYR/XIO_catalogues/output_hard_4.txt"
data = np.genfromtxt(file)
ID = np.array(data[:,0])
ra = np.array(data[:,1])
dec = np.array(data[:,2])
net_cts = np.array(data[:,3])
#flux = np.array(data[:,6])
#flux_err = np.array(data[:,8])

print("\nCalculating new errors...")

offaxis = np.zeros(len(ID))
rr90 = np.zeros(len(ID))
err = np.zeros(len(ID))
ra_err = np.zeros(len(ID))
dec_err = np.zeros(len(ID))

for i in range(len(ID)):
	delta_ra = abs(ra_center - ra[i])
	delta_dec = abs(dec_center - dec[i])
	offaxis[i] = lyr.quadratic_sum(delta_ra, delta_dec)
	rr90[i] = PSF_radius(offaxis[i])
	err[i] = Xerr(rr90[i], net_cts[i])
	ra_err[i] = err[i]/np.sqrt(2)
	dec_err[i] = err[i]/np.sqrt(2)


width=12
name = "Xerr_output_hard_4.txt"
LYRpath = "/Users/Alessandro/Documents/post_TESI/codes/LYR/new_err/"
if not os.path.isdir(LYRpath):
	os.mkdir(LYRpath)

out_file = open(LYRpath + name, "w")
out_file.write("#" + '{0: <{width}.8s}'.format("ID", width=width)
				+ '{0: <{width}.8s}'.format("RA", width=width)
				+ '{0: <{width}.8s}'.format("DEC", width=width)
				+ '{0: <{width}.8s}'.format("err", width=width)
				+ '{0: <{width}.8s}'.format("RA_err", width=width)
				+ '{0: <{width}.8s}'.format("DEC_err", width=width)
				+ '{0: <{width}.8s}'.format("net_cts", width=width)
				+ '{0: <{width}.8s}'.format("off-axis", width=width)
				#+ '{0: <{width}.8s}'.format("flux", width=width)
				#+ '{0: <{width}.8s}'.format("flux_err", width=width)
				+ "\n")
out_file.close()

out_file = open(LYRpath + name, "a")
for d in range(len(ID)):
	out_file.write('{0: <{width}.0f}'.format(ID[d], width=width+1)
				   + '{0: <{width}.6f}'.format(ra[d], width=width)
				   + '{0: <{width}.6f}'.format(dec[d], width=width)
				   + '{0: <{width}.3e}'.format(err[d], width=width)
				   + '{0: <{width}.3e}'.format(ra_err[d], width=width)
				   + '{0: <{width}.3e}'.format(dec_err[d], width=width)
				   + '{0: <{width}.2f}'.format(net_cts[d], width=width)
				   + '{0: <{width}.2f}'.format(offaxis[d]*60, width=width)
				   #+ '{0: <{width}.3e}'.format(flux[d], width=width)
				   #+ '{0: <{width}.3e}'.format(flux_err[d], width=width)
				   + "\n")
out_file.close()

print("... done.\n")


t_ = np.arange(0,15,0.2)
t = (t_*60)/3600
rpsf=[]
for i in range(len(t)):
	rpsf_tmp = PSF_radius(t[i])
	rpsf = np.append(rpsf,rpsf_tmp)


fig,ax=plt.subplots()

ax.plot(t_,rpsf*3600, color='mediumblue', zorder=3, label='Hickox+06')
ax.plot(offaxis*60, rr90*3600, color='darkorange', alpha=0.8, fillstyle='none', lw=0,
		marker='o', zorder=4, label='Data')

ax.set_xlabel('offaxis [arcmin]')
ax.set_ylabel('PSF radius [arcsec]')
#ax.set_yscale('log')
ax.set_title('EER at 90% (E=1.5 keV)')
ax.grid(ls=':', color='grey', alpha=0.5, zorder=0)
ax.legend()

#fig.savefig('xray_errors_puc09.png', bbox_inches="tight", dpi=250)

plt.show()
