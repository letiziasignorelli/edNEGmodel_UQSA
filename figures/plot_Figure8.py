import uncertainpy as un
import numpy as np
import matplotlib.pyplot as plt
from prettyplot import *
import scipy.integrate as integrate
import os

### FIGURE8 ###
data = un.Data()
data.load(filename='../data/simulation_outputs/UQSA_pathological/K_se_5/K_se.h5')
K_se_mean = data['K_se']['mean']
K_se_p5 = data['K_se']['percentile_5']
K_se_p95 = data['K_se']['percentile_95']
K_se_sd = np.sqrt(data['K_se']['variance'])
K_se_var = data['K_se']['variance']
K_se_t = data['K_se']['time']*1e-3

params = [r'$g_\mathrm{Na}$', r'$g_\mathrm{DR}$', r'$g_\mathrm{Ca}$', r'$g_\mathrm{AHP}$', r'$g_\mathrm{C}$']

set_style("seaborn-paper") 
fig = plt.figure()
axesA = plt.subplot(221)
axesB = plt.subplot(222, sharex=axesA)
axesC = plt.subplot(223, sharex=axesA)
axesD = plt.subplot(224, sharex=axesA, sharey=axesB)

fig.set_figwidth(8)
fig.set_figheight(10)

# Panel A
ax = prettyPlot(K_se_t, K_se_mean,
                        color='#8c564b',
                        nr_colors=7,
                        palette="tab10",
                        ax=axesA)

colors = get_current_colormap()
ax.fill_between(K_se_t,
                    K_se_p5,
                    K_se_p95,
                    color=colors[6],
                    linewidth=0)

# Panel B
for i in range(len(params)):
    prettyPlot(K_se_t[1:], data['K_se']['sobol_total'][i,1:],
                            color=i,
                            nr_colors=5,
                            palette="tab10",
                            ax=axesB)

# Panel C
for i in range(len(params)):
    weighted_sobol = data['K_se']['sobol_total'][i,:] * K_se_sd
    prettyPlot(K_se_t[:],weighted_sobol,
                            color=i,
                            nr_colors=5,
                            palette="tab10",
                            ax=axesC)

# Panel D
for i in range(len(params)):
    num = data['K_se']['sobol_total'][i,:] * K_se_var
    general_sobol = []
    for k in range(2,len(num)):
        general_sobol.append(integrate.simpson(num[:k], K_se_t[:k]) / integrate.simpson(K_se_var[:k], K_se_t[:k]))
    # print(general_sobol[-1])
    prettyPlot(K_se_t[2:],general_sobol,
                            color=i,
                            nr_colors=5,
                            palette="tab10",
                            ax=axesD)
    
# Set legend
axesC.legend(labels=params, bbox_to_anchor=(2.02,1.27), fontsize=13, ncol=6, frameon=False, columnspacing=3.2)
axesA.legend([r'$\mathbb{E}[\mathrm{[K^+]_{se}}]$',r'$I_{0.9}$'], fontsize=13, frameon=False, bbox_to_anchor=(1,0.3))

# Axes limits
axesA.set_xlim([min(K_se_t), max(K_se_t)])
axesB.set_xlim([min(K_se_t), max(K_se_t)])
axesC.set_xlim([min(K_se_t), max(K_se_t)])
axesD.set_xlim([min(K_se_t), max(K_se_t)])
axesD.set_ylim([0,1])

# Axes labels
axesA.set_ylabel(r'$\mathrm{[K^+]_{se}}$ $\mathrm{[mM]}$', size = 13)
axesB.set_ylabel(r'$S_{T_i}$', size = 13)
axesC.set_ylabel(r'$S_{T_i}^W$', size = 13)
axesD.set_ylabel(r'$\mathfrak{S}_{T_i}$', size = 13)

axesC.set_xlabel(r'$\mathrm{time}$ $\mathrm{[s]}$', size = 13)
axesD.set_xlabel(r'$\mathrm{time}$ $\mathrm{[s]}$', size = 13)

# ABC
axesA.text(-0.05, 1.1, 'A', transform=axesA.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesB.text(-0.05, 1.1, 'B', transform=axesB.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesC.text(-0.05, 1.1, 'C', transform=axesC.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesD.text(-0.05, 1.1, 'D', transform=axesD.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')

# Set subtitles
axesA.set_title(r"UQ of $\mathrm{[K^+]_{se}}$", fontsize=15, x=0.5, y=1)
axesB.set_title(r"Total order Sobol’ indices", fontsize=15, x=0.5, y=1)
axesC.set_title(r"Weighted Sobol’ indices", fontsize=15, x=0.5, y=1)
axesD.set_title(r"Generalized Sobol’ indices", fontsize=15, x=0.5, y=1)

plt.tight_layout()
plt.subplots_adjust(wspace=0.25, hspace=0.4)

# Save path
# checking if the directory exist and create it if it doesn't
if not os.path.exists('full_figures'):
    os.makedirs('full_figures')
plt.savefig('full_figures/Figure8.pdf', dpi=600)
# plt.show()
