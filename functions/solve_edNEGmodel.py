from scipy.integrate import solve_ivp
import numpy as np
import time
from edNEGmodel.edNEGmodel import *
from data.initial_values.initial_values import *
import h5py
import os

def solve_edNEGmodel_reference(I_stim, alpha, t_span, stim_start, stim_end, path, solver, t_step, jacobian=False):
        """
        Solves the edNEG model with rescaled units (Sætra et al. 2021, v2.0.0) using the solve_ivp function from scipy.

        Arguments:
            I_stim (float): stimulus current [uA]
            alpha (float): coupling strength
            t_span (tuple(float,float)): duration of simulation [ms]
            stim_start (float): time of stimulus onset [ms]
            stim_end (float): time of stimulus offset [ms]
            path (str): The path to store results
            solver (str): The solver name used to numerically solve the differential equations by solve_ivp
            t_step (float): The maximum time step used for numerical integration during the simulation.
            jacobian (bool, optional, default=False): Flag indicating whether to use the Jacobian matrix for the solver.
        """
        hf = h5py.File(path, 'w')

        # Solve
        if jacobian is False:
            start_time = time.time()
            sol = solve_ivp(dkdt, t_span, k0, method=solver, max_step=t_step, args=(alpha, stim_start, stim_end, I_stim))        
            elapsed_time = round(time.time() - start_time, 5)
            print('Solver = ' + solver + ', t_step = ' + str(t_step))
        else:
            start_time = time.time()
            sol = solve_ivp(dkdt, t_span, k0, method=solver, max_step=t_step, args=(alpha, stim_start, stim_end, I_stim), jac=model_jacobian)        
            elapsed_time = round(time.time() - start_time, 5)
            print('Solver = ' + solver + 'j, t_step = ' + str(t_step))

        print('elapsed time: ', elapsed_time, 'seconds')
        print("-------------------------------------------------")
        print(sol)

        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg  = sol.y
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg)

        phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg = my_cell.membrane_potentials()
        phi_de = np.zeros(len(phi_dg))
        membrane_potentials = [phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg]

        E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn = my_cell.reversal_potentials()
        reversal_potentials = [E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn]

        # Save values in .h5 file
        hf.create_dataset('elapsed_time', data = elapsed_time)
        hf.create_dataset('time' , data = sol.t)
        hf.create_dataset('sol' , data = sol.y)
        hf.create_dataset('phi', data =  membrane_potentials)
        hf.create_dataset('E', data = reversal_potentials)
        hf.create_dataset('nfev', data = sol.nfev)
        hf.create_dataset('njev', data = sol.njev)
        hf.create_dataset('nlu' , data = sol.nlu)


def solve_edNEGmodel_convergence(I_stim, alpha, t_span, stim_start, stim_end, path, solvers, t_steps, jacobian=False):
    """
        Perform convergence analysis for the edNEG model with rescaled units (Sætra et al. 2021, v2.0.0) using the solve_ivp function from scipy.

        Arguments:
            I_stim (float): stimulus current [uA]
            alpha (float): coupling strength
            t_span (tuple(float,float)): duration of simulation [ms]
            stim_start (float): time of stimulus onset [ms]
            stim_end (float): time of stimulus offset [ms]
            path (str): The path to store results
            solvers (list of str): A list og the solver names used to numerically solve the differential equations by solve_ivp
            t_steps (list of float): A list of the maximum time steps used for numerical integration during the simulation.
            jacobian (bool, optional, default=False): Flag indicating whether to use the Jacobian matrix for the solver.
        """
    for solver in solvers:   
        # Save path
        if jacobian is False:
            path = os.path.join(path, solver + '.h5')
            hf = h5py.File(path, 'w') 
        else:
            path = os.path.join(path, solver + 'j.h5')
            hf = h5py.File(path, 'w') 


        for t_step in t_steps:
            # Save group 
            group = hf.create_group(str(t_step)) 
            
            # Solve
            if jacobian is False:
                start_time = time.time()
                sol = solve_ivp(dkdt, t_span, k0, method=solver, max_step=t_step, args=(alpha, stim_start, stim_end, I_stim))        
                elapsed_time = round(time.time() - start_time, 5)
                print('Solver = ' + solver + ', t_step = ' + str(t_step))
            else:
                start_time = time.time()
                sol = solve_ivp(dkdt, t_span, k0, method=solver, max_step=t_step, args=(alpha, stim_start, stim_end, I_stim), jac=model_jacobian)        
                elapsed_time = round(time.time() - start_time, 5)
                print('Solver = ' + solver + 'j, t_step = ' + str(t_step))

            print('elapsed time: ', elapsed_time, 'seconds')
            print("-------------------------------------------------")

            Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg  = sol.y   #(sol[:,i] for i in range(34))       #sol.y 
            my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg)

            phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg = my_cell.membrane_potentials()
            phi_de = np.zeros(len(phi_dg))
            membrane_potentials = [phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg]

            E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn = my_cell.reversal_potentials()
            reversal_potentials = [E_Na_sn, E_Na_sg, E_Na_dn, E_Na_dg, E_K_sn, E_K_sg, E_K_dn, E_K_dg, E_Cl_sn, E_Cl_sg, E_Cl_dn, E_Cl_dg, E_Ca_sn, E_Ca_dn]

            # Save values in .h5 file
            group.create_dataset('elapsed_time', data = elapsed_time)
            group.create_dataset('time' , data = sol.t)
            group.create_dataset('sol' , data = sol.y)
            group.create_dataset('phi', data =  membrane_potentials)
            group.create_dataset('E', data = reversal_potentials)
            group.create_dataset('nfev', data = sol.nfev)
            group.create_dataset('njev', data = sol.njev)
            group.create_dataset('nlu' , data = sol.nlu)

# define differential equations
def dkdt(t,k, alpha, stim_start, stim_end, I_stim):
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg = k
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg)
        
        dNadt_sn, dNadt_se, dNadt_sg, dNadt_dn, dNadt_de, dNadt_dg, dKdt_sn, dKdt_se, dKdt_sg, dKdt_dn, dKdt_de, dKdt_dg, dCldt_sn, dCldt_se, dCldt_sg, dCldt_dn, dCldt_de, dCldt_dg, dCadt_sn, dCadt_se, dCadt_dn, dCadt_de = my_cell.dkdt()
        dndt, dhdt, dsdt, dcdt, dqdt, dzdt = my_cell.dmdt()
        dVsidt, dVsedt, dVsgdt, dVdidt, dVdedt, dVdgdt = my_cell.dVdt()
        
        if t > stim_start and t < stim_end:
            dKdt_sn += I_stim / my_cell.F
            dKdt_se -= I_stim / my_cell.F
        
        return np.array([dNadt_sn, dNadt_se, dNadt_sg, dNadt_dn, dNadt_de, dNadt_dg, dKdt_sn, dKdt_se, dKdt_sg, dKdt_dn, dKdt_de, dKdt_dg, \
            dCldt_sn, dCldt_se, dCldt_sg, dCldt_dn, dCldt_de, dCldt_dg, dCadt_sn, dCadt_se, dCadt_dn, dCadt_de, \
            dndt, dhdt, dsdt, dcdt, dqdt, dzdt, dVsidt, dVsedt, dVsgdt, dVdidt, dVdedt, dVdgdt])

# define analytical jacobian function
def model_jacobian(t,k, alpha, stim_start, stim_end, I_stim):
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg = k
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg)

        return my_cell.edNEG_jacobian(dense=True)