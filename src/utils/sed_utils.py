import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

class SEDUtils:

    def compute_luminosity(flux, redshift):
        cosmo = FlatLambdaCDM(H0=70, Om0=0.29, Tcmb0=2.725)
        dl = cosmo.luminosity_distance(redshift).value  # Mpc
        lum_dist = dl * 3.086e24  # convert to cm
        L = 4 * np.pi * lum_dist**2 * flux
        return L

    def compute_flux(mAB: float, wavelength: float) -> float:
        JY = 1e-23  # Janskys
        C_LIGHT = 299792458.0  # m/s
        wavelength_m = wavelength * 1e-10
        f_nu = 3631 * JY * 10 ** (-0.4 * mAB)
        nu = C_LIGHT / wavelength_m
        nuf_nu = nu * f_nu
        return nuf_nu
    
    
    
