import numpy as np
import matplotlib.pyplot as plt
import os
import time

import LYR_functions as lyr
from LYR_input_parameters import *

toc = time.time()

if save_images == True:
	if not os.path.isdir(path_images):
		os.mkdir(path_images)

if not os.path.isdir(path_output):
	os.mkdir(path_output)

from LYR_LR_compl_part4 import *

if save_output == True:
	from LYR_LRRe_plot import *


tic = time.time()
tt = tic - toc
tt_str = 'seconds.'
if tt > 60:
	tt = tt/60
	tt_str = 'minutes.'
elif tt > 3600:
	tt = tt/3600
	tt_str = 'hours.'
	
print('\nComputational time:', round(tt,1), tt_str, '\n')

if save_output == True:
	outinfo = open(path_output+filename_outinfo,"a")
	outinfo.write('\nComputational time: ' + str(round(tt,1)) + ' ' + tt_str + '\n')
	outinfo.close()


