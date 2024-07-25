import sys
sys.path.append("/Users/admin/Documents/GitHub/Matched-Catalog-Analysis/src/utils")

from config_loader import ConfigLoader
config = ConfigLoader("/Users/admin/Documents/GitHub/Matched-Catalog-Analysis/src/configs/config.yaml")

############################
### LYR input parameters ###
############################

import numpy as np
import pandas as pd
print('\n')
print('Reading catalogues and parameters..')


# f(r) searching radius:
r_fr = config.MATCH_SEARCH_RAD # in arcsec
# total(m) searching radius:
r_in_tm = config.MATCH_SEARCH_RAD # in arcsec
# LR searching radius for FR
r_lr = config.MATCH_SEARCH_RAD #in arcsec

# n(m) annulus searching annulus:
r_min_nm = config.BACKGROUND_INNER_RAD   # in arcsec
r_max_nm = config.BACKGROUND_OUTER_RAD

# distributions bins:
distrib_bins = 21

# results:
path_LR = './src/data_lyr/'
path_output = path_LR + 'output/'
path_images = path_LR + 'images/'

# output catalogue errors moltiplicative factor
plus_erf = 1 # wavdetect errors factor, put 1 if there isn't

# sigma finder ('all' == everwy single source; '1sigma' == 1 sigma mean errror):
sigma_out_finder = '1sigma'

# input catalog type(0 = static pos. error (put 0 below is there are no errors); 1 = usual errors;  2 = two catalogue with no errors):
input_cat_type = 0

# catalogues parameters:
noxy = -99. # coords of undetected sources
nomag = -99. # mag/fluxes of undetected sources

# use of fluxes of IO catalogue:
delta_f2 = False
delta_f2_fluxerr = 0 # 0: static error; 1: every single error for the input cat

# python interactions
iwp = False
cursors = False
save_images = False
save_output = True

# sdding string to all output images
add_str = 'Final_XJ'
# adding string to images titles
add_title = add_str

# out catalog and info
filename_LR = 'LR_' + add_str + '.txt'
filename_outinfo = 'LR_outinfo_' + add_str + '.txt'


############################################################### import data ########
# coordinates in degrees

######################################################################################################## ONIR - X:
# Input catalogue:
file_input = "/Users/admin/Documents/GitHub/Matched-Catalog-Analysis/src/data_lyr/input_catalogs/h20_LYR.csv"

data_input = pd.read_csv(file_input).dropna().to_numpy()
ID_input = np.array(data_input[:,0])
ra_input = np.array(data_input[:,1])
dec_input = np.array(data_input[:,2])
mag_input = np.array(data_input[:,5])


# Output catalogue:
file_output = "/Users/admin/Documents/GitHub/Matched-Catalog-Analysis/src/data_lyr/input_catalogs/cdfs_LYR.csv"


data_output = pd.read_csv(file_output).dropna().to_numpy()
ID_output = np.array(data_output[:,0])
ra_output = np.array(data_output[:,1])
dec_output = np.array(data_output[:,2])
ra_err_output = np.array(data_output[:,3]) / 3600
dec_err_output = np.array(data_output[:,4]) / 3600


#---------------------------------------------------------------------- for sigma selection ------

# Sigma input (ONIR):
data1 = pd.read_csv(file_input).to_numpy()
ID1_input_s = np.array(data1[:,0])
ra1_input_s = np.array(data1[:,1])
dec1_input_s = np.array(data1[:,2])
mag1_input_s = np.array(data1[:,5])

if input_cat_type == 0:
	sigma_input = 0 # 0.3 per cosmos

elif input_cat_type == 1:
	ra1_input_s_err = np.array(data1[:,30])
	dec1_input_s_err = np.array(data1[:,30])

elif input_cat_type == 2:
	mag2_input_s = np.array(data1[:,30])
	ID2_input_s = np.array(data1[:,34])
	ra2_input_s = np.array(data1[:,35])
	dec2_input_s = np.array(data1[:,36])


# Sigma output (X-ray):
ID_output_s = ID_output
ra_output_s = ra_output
dec_output_s = dec_output
ra_err_output_s = ra_err_output
dec_err_output_s = dec_err_output



print('... done.\n')
