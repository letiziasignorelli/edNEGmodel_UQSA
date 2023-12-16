import numpy as np
import uncertainpy as un
from scipy.integrate import solve_ivp
from edNEGmodel.edNEGmodel_params import *
from functions.initial_values import *

def model_setup(t,k, g_Na, g_DR, g_Ca, g_AHP, g_C):
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg = k
        params = np.array([g_Na, g_DR, g_Ca, g_AHP, g_C])
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg, params)
        dNadt_sn, dNadt_se, dNadt_sg, dNadt_dn, dNadt_de, dNadt_dg, dKdt_sn, dKdt_se, dKdt_sg, dKdt_dn, dKdt_de, dKdt_dg, dCldt_sn, dCldt_se, dCldt_sg, dCldt_dn, dCldt_de, dCldt_dg, dCadt_sn, dCadt_se, dCadt_dn, dCadt_de = my_cell.dkdt()        
        dndt, dhdt, dsdt, dcdt, dqdt, dzdt = my_cell.dmdt()        
        dVsidt, dVsedt, dVsgdt, dVdidt, dVdedt, dVdgdt = my_cell.dVdt()  

        if t > stim_start and t < stim_end:
            dKdt_sn += I_stim / my_cell.F
            dKdt_se -= I_stim / my_cell.F

        return np.array([dNadt_sn, dNadt_se, dNadt_sg, dNadt_dn, dNadt_de, dNadt_dg, dKdt_sn, dKdt_se, dKdt_sg, dKdt_dn, dKdt_de, dKdt_dg, \
            dCldt_sn, dCldt_se, dCldt_sg, dCldt_dn, dCldt_de, dCldt_dg, dCadt_sn, dCadt_se, dCadt_dn, dCadt_de, \
            dndt, dhdt, dsdt, dcdt, dqdt, dzdt, dVsidt, dVsedt, dVsgdt, dVdidt, dVdedt, dVdgdt])


def model_jacobian_setup(t,k,g_Na, g_DR, g_Ca, g_AHP, g_C):
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg = k
        params = np.array([g_Na, g_DR, g_Ca, g_AHP, g_C])
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg, params)
        
        return my_cell.edNEG_jacobian(dense=True)

def compute_first_AP(t,variable,interval):
    index_left = np.where(t >= interval[0])[0][0]
    index_right = np.where(t >= interval[1])[0][0]

    local_time = t[index_left:index_right]
    local_variable = variable[index_left:index_right]
    firing = False
    i = 0
    control_firing = False

    while i in range(len(local_time)) and control_firing == False:
        if firing == False and local_variable[i] > -10:
            firing = True
        if firing==True and control_firing == False and local_variable[i-1]<local_variable[i] and local_variable[i]>local_variable[i+1]:
            control_firing = True
        i +=1

    t_val = local_time[i]
    return t_val

def phi_msn(g_Na = 30.  ,       
            g_DR = 15.   ,        
            g_Ca = 11.8  ,        
            g_AHP = .8   ,      
            g_C = 15. ):  
        
        params = [g_Na, g_DR, g_Ca, g_AHP, g_C]
        sol = solve_ivp(model_setup, t_span, k0, method='Radau', args=params, max_step=1e1, jac=model_jacobian_setup)
        Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg  = sol.y
        my_cell = edNEGmodel(T, Na_sn, Na_se, Na_sg, Na_dn, Na_de, Na_dg, K_sn, K_se, K_sg, K_dn, K_de, K_dg, Cl_sn, Cl_se, Cl_sg, Cl_dn, Cl_de, Cl_dg, Ca_sn, Ca_se, Ca_dn, Ca_de, X_sn, X_se, X_sg, X_dn, X_de, X_dg, alpha, cbK_se, cbK_sg, cbK_de, cbK_dg, cbCa_sn, cbCa_dn, n, h, s, c, q, z, V_sn, V_se, V_sg, V_dn, V_de, V_dg, cM_sn, cM_se, cM_sg, cM_dn, cM_de, cM_dg, params)
        phi_sn, phi_se, phi_sg, phi_dn, phi_de, phi_dg, phi_msn, phi_mdn, phi_msg, phi_mdg = my_cell.membrane_potentials()
        
        time = sol.t
        values = phi_msn

        # START OF DEPOLARIZATION BLOCK
        # Time when the neuron enters the depolarization block.
        # First time when the derivative of the membrane potential reaches a constant function below a tolerance of 1e-1
        var_diff = np.diff(values)
        index_final = []
        for i in range(len(time)-1):
            index_array = []
            while i < len(time)-1 and var_diff[i] < 1e-1:
                index_array.append(i)
                i +=1
            if len(index_array) > len(index_final):
                index_final = index_array
        start_depolarization_block = time[index_final[0]]

        # TIME BEFORE FIRST ACTION POTENTIAL
        time_before_first_spike = compute_first_AP(time,phi_msn,[stim_start,stim_end]) - stim_start

        info = {"stimulus_start": stim_start, 
                "stimulus_end": stim_end, 
                "start_depolarization_block" : start_depolarization_block,
                "time_before_first_spike": time_before_first_spike}

        return time, values, info


'''
Pathological conditions
'''
I_stim = 20e-5 # [uA]
alpha = 2
t_dur = 6e3       # [ms]
t_span = (0, t_dur)
stim_start = 0.2e3
stim_end = 5.5e3 


'''
Run UQSA
'''
if __name__ == '__main__': 

    # Initialize the model
    model = un.Model(run=phi_msn, 
                     labels=["Time (ms)", "Membrane potential (mV)"],
                     ignore=True)
    
    def start_depolarization_block(time, values, info):
        return None, info["start_depolarization_block"]
    
    def time_before_first_spike(time, values, info):
        return None, info["time_before_first_spike"]
    
    feature_list = [start_depolarization_block, time_before_first_spike]

    '''
    sigma = 1%
    '''
    # Define a parameter dictionary
    parameters = {  "g_Na" : 30.,  
                    "g_DR" : 15., 
                    "g_Ca" : 11.8,
                    "g_AHP" : .8,
                    "g_C" : 15.
                    }

    # Create the parameters
    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a 1% interval around their fixed value
    parameters.set_all_distributions(un.uniform(0.02))

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                    parameters=parameters,
                                    features=feature_list)
    
    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10, rosenblatt=True, data_folder = "data/simulation_outputs/UQSA_pathological/phi_msn_1")

    '''
    sigma = 5%
    '''
    # Define a parameter dictionary
    parameters = {  "g_Na" : 30.,  
                    "g_DR" : 15., 
                    "g_Ca" : 11.8,
                    "g_AHP" : .8,
                    "g_C" : 15.
                    }

    # Create the parameters
    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a 5% interval around their fixed value
    parameters.set_all_distributions(un.uniform(0.1))

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                    parameters=parameters,
                                    features=feature_list)
    
    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10, rosenblatt=True, data_folder = "data/simulation_outputs/UQSA_pathological/phi_msn_5")

    '''
    sigma = 10%
    '''
    # Define a parameter dictionary
    parameters = {  "g_Na" : 30.,  
                    "g_DR" : 15., 
                    "g_Ca" : 11.8,
                    "g_AHP" : .8,
                    "g_C" : 15.
                    }

    # Create the parameters
    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a 10% interval around their fixed value
    parameters.set_all_distributions(un.uniform(0.2))

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                    parameters=parameters,
                                    features=feature_list)
    
    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10, rosenblatt=True, data_folder = "data/simulation_outputs/UQSA_pathological/phi_msn_10")