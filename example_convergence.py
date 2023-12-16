import numpy as np
from edNEGmodel.edNEGmodel import *
from functions.initial_values import *
from functions.solve_edNEGmodel import *
import os
      
print('*************************************************')
print('*        Convergence analysis starting          *')
print('*************************************************')

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
# checking if the directory exist and create it if it doesn't
savepath = "data/simulation_outputs/convergence_physiological"
if not os.path.exists(savepath):
    os.makedirs(savepath)

# Solve Explicit solvers
solvers = ['RK23', 'RK45', 'DOP853']            
t_steps= [0.1, 0.05, 0.025, 0.0125]   
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=False)

# Solve implicit solvers with analytical Jacobian
solvers = ['LSODA', 'BDF']
t_steps= [np.inf, 1e2, 1e1, 1, 0.1, 0.05, 0.025, 0.0125] 
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=True)

solvers = ['Radau']    
t_steps= [np.inf, 1e2, 1e1, 1, 0.1]  
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=True)

# Solve implicit solvers with FD approximated Jacobian
solvers = ['LSODA', 'BDF']
t_steps= [np.inf, 1e2, 1e1, 1, 0.1, 0.05, 0.025, 0.0125] 
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=False)

solvers = ['Radau']    
t_steps= [np.inf, 1e2, 1e1, 1, 0.1]      
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=False)


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
# checking if the directory exist and create it if it doesn't
savepath = "data/simulation_outputs/convergence_pathological"
if not os.path.exists(savepath):
    os.makedirs(savepath)

# Solve Explicit solvers
solvers = ['RK23', 'RK45', 'DOP853']            
t_steps= [0.1, 0.05, 0.025, 0.0125]   
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=False)

# Solve implicit solvers with analytical Jacobian
solvers = ['LSODA', 'BDF']
t_steps= [np.inf, 1e2, 1e1, 1, 0.1, 0.05, 0.025, 0.0125] 
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=True)

solvers = ['Radau']    
t_steps= [np.inf, 1e2, 1e1, 1, 0.1]  
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=True)

# Solve implicit solvers with FD approximated Jacobian
solvers = ['LSODA', 'BDF']
t_steps= [np.inf, 1e2, 1e1, 1, 0.1, 0.05, 0.025, 0.0125] 
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=False)

solvers = ['Radau']    
t_steps= [np.inf, 1e2, 1e1, 1, 0.1]      
solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, savepath, solvers, t_steps, jacobian=False)