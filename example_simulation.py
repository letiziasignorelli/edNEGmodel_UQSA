import numpy as np
from edNEGmodel.edNEGmodel import *
from data.initial_values.initial_values import *
from functions.solve_edNEGmodel import *
import os

print('*************************************************')
print('*        Computing referece simulation          *')
print('*************************************************')

'''
Resting state (no stimulus current)
'''
I_stim = 0.             # [uA]
alpha = 2
t_dur = 1e4             # [ms]
t_span = (0, t_dur)
stim_start = 0.
stim_end = 0.

# Save path
# checking if the directory exist and create it if it doesn't
savepath = "data/simulation_outputs/example"
file_name = "restingstate.h5"
if not os.path.exists(savepath):
    os.makedirs(savepath)
path = os.path.join(savepath, file_name)

# Choose simulation setting
solver = 'Radau'
t_step = 1e1

# Solve
solve_edNEGmodel_reference(I_stim, alpha, t_span, stim_start, stim_end, path, solver, t_step, jacobian=True)


'''
Physiological conditions
'''
I_stim = 7.2e-5         # [uA]
alpha = 2
t_dur = 1e4             # [ms]
t_span = (0, t_dur)
stim_start = 1e3
stim_end = 5e3         

# Save path
file_name = "physiological.h5"
path = os.path.join(savepath, file_name)
# Solve
solve_edNEGmodel_reference(I_stim, alpha, t_span, stim_start, stim_end, path, solver, t_step, jacobian=True)


'''
Pathological conditions
'''
I_stim = 15e-5          # [uA]
alpha = 2
t_dur = 1e4             # [ms]
t_span = (0, t_dur)
stim_start = 1e3
stim_end = 8e3  

# Save path
file_name = "pathological.h5"
path = os.path.join(savepath, file_name)
# Solve
solve_edNEGmodel_reference(I_stim, alpha, t_span, stim_start, stim_end, path, solver, t_step, jacobian=True)