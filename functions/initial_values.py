import numpy as np
import pkg_resources

filename = pkg_resources.resource_filename('data', 'initial_values/initial_values.npz')
data = np.load(filename)

T = 309.14 #[K]

# initial membrane potentials [mV]
phi_msn0 = data['phi_msn'] * 1e3
phi_msg0 = data['phi_msg'] * 1e3
phi_mdn0 = data['phi_mdn'] * 1e3
phi_mdg0 = data['phi_mdg'] * 1e3

# initial volumes [cm**3]
V_sn0 = data['V_sn'] * 1e6
V_se0 = data['V_se'] * 1e6
V_sg0 = data['V_sg'] * 1e6
V_dn0 = data['V_dn'] * 1e6
V_de0 = data['V_de'] * 1e6
V_dg0 = data['V_dg'] * 1e6

# initial amounts of ions [nmol]
Na_sn0 = data['Na_sn'] * 1e9
Na_se0 = data['Na_se'] * 1e9
Na_sg0 = data['Na_sg'] * 1e9
K_sn0 = data['K_sn'] * 1e9
K_se0 = data['K_se'] * 1e9
K_sg0 = data['K_sg'] * 1e9
Cl_sn0 = data['Cl_sn'] * 1e9
Cl_se0 = data['Cl_se'] * 1e9
Cl_sg0 = data['Cl_sg'] * 1e9
Ca_sn0 = data['Ca_sn'] * 1e9
Ca_se0 = data['Ca_se'] * 1e9

Na_dn0 = data['Na_dn'] * 1e9
Na_de0 = data['Na_de'] * 1e9
Na_dg0 = data['Na_dg'] * 1e9
K_dn0 = data['K_dn']  * 1e9
K_de0 = data['K_de'] * 1e9
K_dg0 = data['K_dg'] * 1e9
Cl_dn0 = data['Cl_dn'] * 1e9
Cl_de0 = data['Cl_de'] * 1e9
Cl_dg0 = data['Cl_dg'] * 1e9
Ca_dn0 = data['Ca_dn'] * 1e9
Ca_de0 = data['Ca_de'] * 1e9

# intial gating variables
n0 = data['n']
h0 = data['h']
s0 = data['s']
c0 = data['c']
q0 = data['q']
z0 = data['z']

# baseline ion concentrations [nmol/cm**3 ]
cbK_se = 3.082 * 1e3
cbK_sg = 99.959 * 1e3
cbK_de = 3.082 * 1e3
cbK_dg = 99.959 * 1e3
cbCa_sn = 0.01 * 1e3
cbCa_dn = 0.01 * 1e3

# residual charges [nmol]
res_sn = phi_msn0*3e-2*616e-12/9.648e4 *1e6
res_sg = phi_msg0*3e-2*616e-12/9.648e4  *1e6
res_se = res_sn+res_sg
res_dn = phi_mdn0*3e-2*616e-12/9.648e4  *1e6
res_dg = phi_mdg0*3e-2*616e-12/9.648e4  *1e6
res_de = res_dn+res_dg

X_sn = Na_sn0 + K_sn0 - Cl_sn0 + 2*Ca_sn0 - res_sn
X_se = Na_se0 + K_se0 - Cl_se0 + 2*Ca_se0 + res_se
X_sg = Na_sg0 + K_sg0 - Cl_sg0 - res_sg
X_dn = Na_dn0 + K_dn0 - Cl_dn0 + 2*Ca_dn0 - res_dn
X_de = Na_de0 + K_de0 - Cl_de0 + 2*Ca_de0 + res_de
X_dg = Na_dg0 + K_dg0 - Cl_dg0 - res_dg

# residual mass [nmol/cm**3]
cM_sn = (Na_sn0 + K_sn0 + Cl_sn0 + Ca_sn0)/V_sn0
cM_se = (Na_se0 + K_se0 + Cl_se0 + Ca_se0)/V_se0 
cM_sg = (Na_sg0 + K_sg0 + Cl_sg0)/V_sg0
cM_dn = (Na_dn0 + K_dn0 + Cl_dn0 + Ca_dn0)/V_dn0
cM_de = (Na_de0 + K_de0 + Cl_de0 + Ca_de0)/V_de0 
cM_dg = (Na_dg0 + K_dg0 + Cl_dg0)/V_dg0 

k0 = [Na_sn0, Na_se0, Na_sg0, Na_dn0, Na_de0, Na_dg0, K_sn0, K_se0, K_sg0, K_dn0, K_de0, K_dg0, Cl_sn0, Cl_se0, Cl_sg0, Cl_dn0, Cl_de0, Cl_dg0, Ca_sn0, Ca_se0, Ca_dn0, Ca_de0, n0, h0, s0, c0, q0, z0, V_sn0, V_se0, V_sg0, V_dn0, V_de0, V_dg0]

