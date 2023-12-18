import matplotlib.pyplot as plt
from prettyplot import *
import h5py
import os
import numpy as np

# TYPE CHOICE
'''
type choice =   'sol'   : solution (Number of ions, Gating variables, Volumes)
                'phi'   : membrane potentials
                'E'     : reversal potentials
'''

# VAR CHOICE
'''
for type = sol
    --> choose: 'Na_sn', 'Na_se', 'Na_sg', 'Na_dn', 'Na_de', 'Na_dg', 'K_sn', 'K_se', 'K_sg', 'K_dn', 'K_de', 'K_dg', 'Cl_sn', 'Cl_se', 
                'Cl_sg', 'Cl_dn', 'Cl_de', 'Cl_dg', 'Ca_sn', 'Ca_se', 'Ca_dn', 'Ca_de', 'n', 'h', 's', 'c', 'q', 'z', 'V_sn', 'V_se', 'V_sg', 'V_dn', 'V_de', 'V_dg'
for type = phi
    --> choose: 'phi_sn', 'phi_se', 'phi_sg', 'phi_dn', 'phi_de', 'phi_dg', 'phi_msn', 'phi_mdn', 'phi_msg', 'phi_mdg'
for type = E
    --> choose: 'E_Na_sn', 'E_Na_sg', 'E_Na_dn', 'E_Na_dg', 'E_K_sn', 'E_K_sg', 'E_K_dn', 'E_K_dg', 'E_Cl_sn', 'E_Cl_sg', 'E_Cl_dn', 'E_Cl_dg', 'E_Ca_sn', 'E_Ca_dn'
'''


### FIGURE2 ###
types = ['phi', 'sol', 'sol']
vars = ['phi_msn', 'K_se', 'V_se']
sim_folder = '../data/simulation_outputs/'
solvers = ['UQSA_physiological/nominal','UQSA_pathological/nominal']

set_style("seaborn-paper")
fig = plt.figure()
axesA = plt.subplot(321)
axesB = plt.subplot(322, sharey=axesA)
axesC = plt.subplot(323)
axesD = plt.subplot(324)
axesE = plt.subplot(325)
axesF = plt.subplot(326, sharey=axesE)
axes = np.asarray([[axesA, axesB], [axesC, axesD], [axesE, axesF]])

fig.set_figwidth(10)
fig.set_figheight(7.5)

for i in range(0, 2*len(vars)):
    nx = i % 2
    ny = int(np.floor(i/float(2)))

    if i < 2*len(vars):
        path = os.path.join(sim_folder, solvers[nx] + '.h5')
        hf = h5py.File(path, 'r') 
        t_ref = hf['time'][()]*1e-3

        sol_index = switch(types[ny],vars[ny])  
        variable_ref = hf[types[ny]][(sol_index)]  

        if vars[ny] in concentrations:
            sol_index1 = switch('sol','V_se')  
            variable_ref = variable_ref / hf['sol'][(sol_index1)] *1e-3

        if vars[ny] in volumes:
            variable_ref = (variable_ref - variable_ref[0])/variable_ref[0]*100

        prettyPlot(t_ref, variable_ref,
                    color=ny,
                    nr_colors=2*len(vars),
                    ax=axes[ny][nx],
                    palette="tab10")

        axes[ny][nx].set_xlim([min(t_ref), max(t_ref)])

# Axes limits
axesA.set_xlim([0,6])
axesD.set_xlim([0,6])

# Axes labels
axesA.set_ylabel(r'$\phi_\mathrm{msn}$ $\mathrm{[mV]}$', size = 13)
axesC.set_ylabel(r'$\mathrm{[K^+]_{se}}$ $\mathrm{[mM]}$', size = 13)
axesE.set_ylabel(r'$\Delta \mathrm{V_{se}}$ $[\%]$', size = 13)

axesE.set_xlabel(r'$\mathrm{time}$ $\mathrm{[s]}$', size = 13)
axesF.set_xlabel(r'$\mathrm{time}$ $\mathrm{[s]}$', size = 13)

# Set titles
axesA.set_title('Physiological', size=20, pad=12)
axesB.set_title('Pathological', size=20, pad=12)

plt.tight_layout()
plt.subplots_adjust(top=0.9, wspace=0.25, hspace=0.45)
plt.savefig('full_figures/Figure2.png', dpi=600)
# plt.show()



