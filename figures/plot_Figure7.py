import uncertainpy as un
import numpy as np
import matplotlib.pyplot as plt
from prettyplot import * 

### FIGURE7 ###
data = un.Data()
data.load(filename='../data/simulation_outputs/UQSA_pathological/phi_msn_5/phi_msn.h5')
sobol_1_sdp = data['start_depolarization_block']['sobol_first']
sobol_tot_sdp = data['start_depolarization_block']['sobol_total']
sobol_1_tbfs = data['time_before_first_spike']['sobol_first']
sobol_tot_tbfs = data['time_before_first_spike']['sobol_total']

xlab = [r'$g_\mathrm{Na}$', r'$g_\mathrm{DR}$', r'$g_\mathrm{Ca}$', r'$g_\mathrm{AHP}$', r'$g_\mathrm{C}$']

set_style("seaborn-paper")
fig = plt.figure()
axesA = plt.subplot(221)
axesB = plt.subplot(222, sharex=axesA, sharey=axesA)
axesC = plt.subplot(223, sharex=axesA, sharey=axesA)
axesD = plt.subplot(224, sharex=axesA, sharey=axesA)

fig.set_figwidth(7)
fig.set_figheight(6)

width = 0.2
index = np.arange(1, len(data.uncertain_parameters)+1)*width

# Panel A - C
prettyBar(sobol_1_sdp,
            palette="tab10",
            xlabels=xlab,
            index=index,
            ax=axesA)
prettyBar(sobol_tot_sdp,
            palette="tab10",
            xlabels=xlab,
            index=index,
            ax=axesC)

# Panel B - D
prettyBar(sobol_1_tbfs,
            palette="tab10",
            xlabels=xlab,
            index=index,
            ax=axesB)
prettyBar(sobol_tot_tbfs,
            palette="tab10",
            xlabels=xlab,
            index=index,
            ax=axesD)

axesA.set_ylim([0, 1])

axesA.set_ylabel(r"$S_i$", size=13)
axesC.set_ylabel(r"$S_{T_i}$", size=13)

axesA.set_title(r"$T_\mathrm{sDP}$", size=20, pad=17)
axesB.set_title(r'$T_\mathrm{bFAP}$', size=20, pad=17)

# ABC
# axesA.text(-0.05, 1.2, 'A', transform=axesA.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
# axesB.text(-0.05, 1.2, 'B', transform=axesB.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
# axesC.text(-0.05, 1.2, 'C', transform=axesC.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
# axesD.text(-0.05, 1.2, 'D', transform=axesD.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
# axesE.text(-0.05, 1.2, 'E', transform=axesE.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
# axesF.text(-0.05, 1.2, 'F', transform=axesF.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')

plt.tight_layout()
plt.subplots_adjust(top=0.9, wspace=0.25, hspace=0.4)
plt.savefig('full_figures/Figure7.png', dpi=600)
# plt.show()



