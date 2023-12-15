T = 309.14 #[K]

# initial membrane potentials [mV]
phi_msn0 = -66.93391011908824
phi_msg0 = -83.90412349979285
phi_mdn0 = -66.93175609621206
phi_mdg0 = 83.89970092819222

# initial volumes [cm**3]
V_sn0 = 1.437e-09
V_se0 = 7.185e-10
V_sg0 = 1.437e-09
V_dn0 = 1.437e-09
V_de0 = 7.185e-10
V_dg0 = 1.437e-09

# initial amounts of ions [nmol]
Na_sn0 = 2.693135361094143e-05
Na_se0 = 0.00010227456844383901
Na_sg0 = 2.0820741489654975e-05
K_sn0 = 0.00019839631683593885
K_se0 = 2.5434848036884403e-06
K_sg0 = 0.00014537832547008566
Cl_sn0 = 1.0267828109597383e-05
Cl_se0 = 9.476324371500056e-05
Cl_sg0 = 8.124814212251761e-06
Ca_sn0 = 1.4370000663073967e-08
Ca_se0 = 7.902454912405598e-07

Na_dn0 = 2.6945093581068462e-05
Na_de0 = 0.00010225667114059411
Na_dg0 = 2.081717174243654e-05
K_dn0 = 0.00019838276831431778
K_de0 = 2.550653957423747e-06
K_dg0 = 0.00014538245061980097
Cl_dn0 = 1.0268131705275308e-05
Cl_de0 = 9.475282229768452e-05
Cl_dg0 = 8.125368767639134e-06
Ca_dn0 = 1.442628054228976e-08
Ca_de0 = 7.903982275609302e-07

# intial gating variables
n0 = 0.00030640329262244734
h0 = 0.9993074990050067
s0 = 0.007662867534221651
c0 = 0.005653177185282168
q0 = 0.01169467703634858
z0 = 1.0

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

