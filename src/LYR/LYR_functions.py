import numpy as np


################################################################ generic fuctions:
def quadratic_sum(x,y):
	return np.sqrt(np.power(x,2) + np.power(y,2))


#annulus area:
def aarea(rmin, rmax):
	return np.pi*(np.power(rmax,2) - np.power(rmin,2))


def acircle(r):
	return np.pi*np.power(r,2)


#from barplot to histogram:
def bartohist(x,y):
	y_len = int(np.sum(np.around(y,decimals=0)))
	new_y = np.array(np.around(y,decimals=0), dtype=int)
	new_x = np.zeros(y_len)
	a=b=j=0
	index=True
	while(index==True):
		if a < new_y[j]:
			new_x[b] = x[j]
			b=b+1
			a=a+1
		else:
			j=j+1
			a=0
		if b == len(new_x):
			index = False
	return new_x


def gaussian(x,mu,sigma,norm):
	return norm*np.exp(-(x-mu)**2/(2*np.power(sigma,2)))



################################################################ specific fuctions:

#f(r,c): angular separation prob. distribution, where
#sigma_o: 1sigma position error in onir srcs; sigma_x: 1sigma position error in X srcs; r=sqrt(ra_off**2+dec_off**2)
def FR(sigma_out, sigma_in, deltaRA, deltaDEC):
	sigma_ = quadratic_sum(sigma_in, sigma_out)
	r = quadratic_sum(deltaRA, deltaDEC)
	fr = (1/(2*np.pi*np.power(sigma_,2)))*np.exp(-(np.power(r,2))/(2*np.power(sigma_,2)))
	return r, fr


#LR: likelihood ratio
def LR(fr,nm,qm):
	return qm*fr/nm

#LR: likelihood ratio + fdf
def LR_(fr,nm,qm,fdf):
	return qm*fr*fdf/nm


#Re: realibiliy
def Re(summLR, LR_single, Q):
	Re = LR_single/(summLR+(1-Q))
	return Re


# f(delta_f): multiplicative function to take into account the output flux
def Delta_flux(flux_out, flux_in, flux_err_out, flux_err_in):
	#delta_flux = abs(flux_out - flux_in)
	#err_ = quadratic_sum(flux_err_out,flux_err_in)
	#flux_out = np.log10(flux_out)
	#flux_in = np.log10(flux_in)

	delta_flux = np.log10(abs(flux_out - flux_in))
	err_ = np.log10(quadratic_sum(flux_err_out,flux_err_in))
	#fdf = (1/(2*np.pi*np.power(err_,2)))*np.exp(-(np.power(delta_flux,2))/(2*np.power(err_,2)))
	fdf = np.exp(-(np.power(delta_flux,2))/(2*np.power(err_,2)))
	return delta_flux, fdf

