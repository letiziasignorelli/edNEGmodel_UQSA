import matplotlib.pyplot as plt
from prettyplot import *
import h5py
import os
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


### TABLE8 ###
savepath = '../data/simulation_outputs/factor_fixing'
file_name = "factor_fixing.h5"
path = os.path.join(savepath, file_name)
hf = h5py.File(path, 'r')

# Open results
output_samples = hf['output_samples'][()]
first_order_indices = hf['S1'][()]
total_order_indices = hf['ST'][()]

output_list = [r'$\phi_\mathrm{msn}$', r'$\phi_\mathrm{mdn}$', r'$\phi_\mathrm{msg}$', r'$\phi_\mathrm{mdg}$', r'$\mathrm{[K^+]_{se}}$', r'$\mathrm{[K^+]_{de}}$']

with open('full_figures/Table3.txt', 'w') as file:
    print(output_list, file=file)
    print(total_order_indices, file=file)