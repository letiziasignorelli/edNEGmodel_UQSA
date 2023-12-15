from edNEGmodel.edNEGmodel import *
from data.initial_values.initial_values import *
from functions.solve_edNEGmodel import *

'''
Physiological conditions
'''
I_stim = 8e-5 # [uA]
alpha = 2
t_dur = 6e3       # [ms]
t_span = (0, t_dur)
stim_start = 0.2e3
stim_end = 5.5e3 

# Save path
# checking if the directory exist and create it if it doesn't
savepath = "data/simulation_outputs/UQSA_physiological"
file_name = "nominal.h5"
if not os.path.exists(savepath):
    os.makedirs(savepath)
path = os.path.join(savepath, file_name)
hf = h5py.File(path, 'w')

# Choose simulation setting
solver = 'Radau'
t_step = 1e3

# Solve
sol = solve_ivp(dkdt, t_span, k0, method=solver, max_step=t_step, first_step=1e-1, args=(alpha, stim_start, stim_end, I_stim), jac=model_jacobian)
Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg  = sol.y
my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg)
phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg = my_cell.membrane_potentials()
phi_de = np.zeros(len(phi_dg))
membrane_potentials = [phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg]

E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn = my_cell.reversal_potentials()
reversal_potentials = [E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn]

# Save values in .h5 file
hf.create_dataset('time' , data = sol.t)
hf.create_dataset('sol' , data = sol.y)
hf.create_dataset('phi', data =  membrane_potentials)
hf.create_dataset('E', data = reversal_potentials)
hf.create_dataset('nfev', data = sol.nfev)
hf.create_dataset('njev', data = sol.njev)
hf.create_dataset('nlu' , data = sol.nlu)


'''
Pathological conditions
'''
I_stim = 20e-5 # [uA]
alpha = 2
t_dur = 6e3       # [ms]
t_span = (0, t_dur)
stim_start = 0.2e3
stim_end = 5.5e3 

# Save path
# checking if the directory exist and create it if it doesn't
savepath = "data/simulation_outputs/UQSA_pathological"
file_name = "nominal.h5"
if not os.path.exists(savepath):
    os.makedirs(savepath)
path = os.path.join(savepath, file_name)
hf = h5py.File(path, 'w')

# Choose simulation setting
solver = 'Radau'
t_step = 1e1

# Solve
start_time = time.time()
sol = solve_ivp(dkdt, t_span, k0, method=solver, max_step=t_step, first_step=1e-1, args=(alpha, stim_start, stim_end, I_stim), jac=model_jacobian)
Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg  = sol.y
my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg)
phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg = my_cell.membrane_potentials()
phi_de = np.zeros(len(phi_dg))
membrane_potentials = [phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg]
elapsed_time = round(time.time() - start_time, 5)
print('elapsed time: ', elapsed_time, 'seconds')
print("-------------------------------------------------")

E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn = my_cell.reversal_potentials()
reversal_potentials = [E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn]

# Save values in .h5 file
hf.create_dataset('time' , data = sol.t)
hf.create_dataset('sol' , data = sol.y)
hf.create_dataset('phi', data =  membrane_potentials)
hf.create_dataset('E', data = reversal_potentials)
hf.create_dataset('nfev', data = sol.nfev)
hf.create_dataset('njev', data = sol.njev)
hf.create_dataset('nlu' , data = sol.nlu)