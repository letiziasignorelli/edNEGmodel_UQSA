from scipy.integrate import solve_ivp
import numpy as np
import os
import h5py
from edNEGmodel.edNEGmodel_params import * 
from data.initial_values.initial_values import *
from SALib.sample.saltelli import sample as ss
from SALib.analyze.sobol import analyze as sa

'''
Table 3: Factor Fixing analysis
'''

# Define the problem setup

# Choose sigma_hat
ic = 0.15

# Setup problem parameters distributions and groups
problem = {
    'num_vars': 16,
    'names': [  'g_Na_leak_n',
                'g_K_leak_n',
                'g_Cl_leak_n',
                'g_Na' ,      
                'g_DR' ,         
                'g_Ca' ,       
                'g_AHP',         
                'g_C'   ,      
                'g_Na_leak_g' ,  
                'g_K_IR' ,    
                'g_Cl_leak_g' , 
                'rho_n',        
                'U_kcc2' ,  
                'U_nkcc1' ,
                'U_Cadec' ,      
                'rho_g'],
    'bounds': [       
                [0.0246 * (1 - ic), 0.0246 * (1 + ic)],
                [0.0245 * (1 - ic), 0.0245 * (1 + ic)],
                [0.1 * (1 - ic), 0.1 * (1 + ic)],
                [30. * (1 - ic), 30. * (1 + ic)],
                [15. * (1 - ic), 15. * (1 + ic)],
                [11.8 * (1 - ic), 11.8 * (1 + ic)], 
                [.8 * (1 - ic), .8 * (1 + ic)],
                [15. * (1 - ic), 15. * (1 + ic)], 
                [.1 * (1 - ic), .1 * (1 + ic)],
                [1.696 * (1 - ic), 1.696 * (1 + ic)],
                [0.05 * (1 - ic), 0.05 * (1 + ic)], 
                [1.87e-4 * (1 - ic), 1.87e-4 * (1 + ic)],
                [1.49e-5 * (1 - ic), 1.49e-5 * (1 + ic)], 
                [2.33e-5 * (1 - ic), 2.33e-5 * (1 + ic)], 
                [0.075 * (1 - ic), 0.075 * (1 + ic)],
                [1.12e-4 * (1 - ic),1.12e-4 * (1 + ic)],
              ],
    'dists':['unif','unif','unif','unif','unif','unif','unif','unif',\
              'unif','unif','unif','unif','unif','unif','unif','unif'],
    'groups': ['leak','leak','leak','dyn','dyn','dyn','dyn','dyn',\
              'leak','leak','leak','leak','leak','leak','leak','leak']
    }

# Define functions to solve the model
def model_setup(t,k,g_Na_leak_n, g_K_leak_n, g_Cl_leak_n, g_Na, g_DR, g_Ca, g_AHP, g_C, g_Na_leak_g, g_K_IR, g_Cl_leak_g, rho_n, U_kcc2, U_nkcc1, U_Cadec, rho_g):
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg = k
        params = np.array([g_Na_leak_n, g_K_leak_n, g_Cl_leak_n, g_Na, g_DR, g_Ca, g_AHP, g_C, g_Na_leak_g, g_K_IR, g_Cl_leak_g, rho_n, U_kcc2, U_nkcc1, U_Cadec, rho_g])
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg, params)
        dNadt_sn, dNadt_se, dNadt_sg, dNadt_dn, dNadt_de, dNadt_dg, dKdt_sn, dKdt_se, dKdt_sg, dKdt_dn, dKdt_de, dKdt_dg, dCldt_sn, dCldt_se, dCldt_sg, dCldt_dn, dCldt_de, dCldt_dg, dCadt_sn, dCadt_se, dCadt_dn, dCadt_de = my_cell.dkdt()        
        dndt, dhdt, dsdt, dcdt, dqdt, dzdt = my_cell.dmdt()        
        dVsidt, dVsedt, dVsgdt, dVdidt, dVdedt, dVdgdt = my_cell.dVdt()  

        return np.array([dNadt_sn, dNadt_se, dNadt_sg, dNadt_dn, dNadt_de, dNadt_dg, dKdt_sn, dKdt_se, dKdt_sg, dKdt_dn, dKdt_de, dKdt_dg, \
            dCldt_sn, dCldt_se, dCldt_sg, dCldt_dn, dCldt_de, dCldt_dg, dCadt_sn, dCadt_se, dCadt_dn, dCadt_de, \
            dndt, dhdt, dsdt, dcdt, dqdt, dzdt, dVsidt, dVsedt, dVsgdt, dVdidt, dVdedt, dVdgdt])

def model_jacobian(t,k,g_Na_leak_n, g_K_leak_n, g_Cl_leak_n, g_Na, g_DR, g_Ca, g_AHP, g_C, g_Na_leak_g, g_K_IR, g_Cl_leak_g, rho_n, U_kcc2, U_nkcc1, U_Cadec, rho_g):
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg = k
        params = np.array([g_Na_leak_n, g_K_leak_n, g_Cl_leak_n, g_Na, g_DR, g_Ca, g_AHP, g_C, g_Na_leak_g, g_K_IR, g_Cl_leak_g, rho_n, U_kcc2, U_nkcc1, U_Cadec, rho_g])
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg, params)
        
        return my_cell.edNEG_jacobian(dense=True)

def model_run(params):     
        sol = solve_ivp(model_setup, t_span, k0, method=solver, args=params, max_step=t_step, first_step=t_step, jac=model_jacobian)
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg  = sol.y[:,-2]
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg, params)
        phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg = my_cell.membrane_potentials()
        return phi_msn, phi_mdn, phi_msg, phi_mdg, K_se, K_de


# Parameters nominal values
# conductances
g_Na_leak_n = 0.0246
g_K_leak_n = 0.0245
g_Cl_leak_n = 0.1   
g_Na = 30.         
g_DR = 15.           
g_Ca = 11.8          
g_AHP = .8         
g_C = 15.          
g_Na_leak_g = .1    
g_K_IR = 1.696     
g_Cl_leak_g =  0.05     

# exchanger strengths
rho_n = 1.87e-4      
U_kcc2 = 1.49e-5    
U_nkcc1 = 2.33e-5 
U_Cadec = 0.075       
rho_g = 1.12e-4  

params = [g_Na_leak_n, g_K_leak_n, g_Cl_leak_n, g_Na, g_DR, g_Ca, g_AHP, g_C, g_Na_leak_g, g_K_IR, g_Cl_leak_g, rho_n, U_kcc2, U_nkcc1, U_Cadec, rho_g]
params_names = ['g_Na_leak_n', 'g_K_leak_n', 'g_Cl_leak_n', 'g_Na', 'g_DR', 'g_Ca', 'g_AHP', 'g_C', 'g_Na_leak_g', 'g_K_IR', 'g_Cl_leak_g', 'rho_n', 'U_kcc2', 'U_nkcc1', 'U_Cadec', 'rho_g']


print('*************************************************')
print('*              Analysis starting                *')
print('*************************************************')

# Rescaled test parameters in resting conditions
I_stim = 0.             # [uA]
alpha = 2
t_dur = 24e4             # [ms]
t_span = (0, t_dur)
stim_start = 0.
stim_end = 0.

# Save path
# checking if the directory exist and create it if it doesn't
savepath = "data/simulation_outputs/factor_fixing"
file_name = "factor_fixing.h5"
if not os.path.exists(savepath):
    os.makedirs(savepath)
path = os.path.join(savepath, file_name)
hf = h5py.File(path, 'w')

# Choose simulation setting
solver = 'Radau'
t_step = 1e1

# Choose number of samples
num_samples = 2**9

# Input parameters sampled froma Uniform +-15% with Saltelli method
param_samples = ss(problem, num_samples, calc_second_order=True)
print('Parameters sample done!')

# Compute output parameters
output_samples= np.array([model_run(inputs) for inputs in param_samples])
print('Outputs sample done!')

# Save results
hf.create_dataset('output_samples', data = output_samples)

# Perform Sobol sensitivity analysis for each output
first_order_indices = []
total_order_indices = []

for i in range(output_samples.shape[1]):     
    sobol_results = sa(problem, output_samples[:,i], calc_second_order=True)
    first_order_indices.append(sobol_results['S1'])
    total_order_indices.append(sobol_results['ST'])
print('Sobol indices computed!')

# Save results
hf.create_dataset('S1', data = first_order_indices)
hf.create_dataset('ST', data = total_order_indices)

