import uncertainpy as un
import numpy as np
import matplotlib.pyplot as plt
from prettyplot import * 
import os
import h5py
from scipy.stats import zscore

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

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel(r'$\hat{\sigma}$ $[\%]$')

### FIGURE6 ###
# Import nominal data
sim_folder = '../data/simulation_outputs/UQSA_pathological'
path = os.path.join(sim_folder, 'nominal.h5')
hf = h5py.File(path, 'r') 
t_ref = hf['time'][()]

sol_index = switch('phi','phi_msn')  
variable_ref = hf['phi'][(sol_index)] 
var_diff = np.diff(variable_ref)
index_final = []
for i in range(len(t_ref)-1):
    index_array = []
    while i < len(t_ref)-1 and var_diff[i] < 1e-1:
        index_array.append(i)
        i +=1
    if len(index_array) > len(index_final):
        index_final = index_array
start_depolarization_block_nominal = np.full(254, t_ref[index_final[0]]*1e-3)
time_before_first_spike_nominal = np.full(254,compute_first_AP(t_ref,variable_ref,[0.2e3,5.5e3]) - 0.2e3)

# Import UQ data
data1 = un.Data()
data1.load(filename='../data/simulation_outputs/UQSA_pathological/phi_msn_1/phi_msn.h5')
data5 = un.Data()
data5.load(filename='../data/simulation_outputs/UQSA_pathological/phi_msn_5/phi_msn.h5')
data10 = un.Data()
data10.load(filename='../data/simulation_outputs/UQSA_pathological/phi_msn_10/phi_msn.h5')

start_depolarization_block1 = data1['start_depolarization_block']['evaluations']*1e-3
start_depolarization_block5 = data5['start_depolarization_block']['evaluations']*1e-3
start_depolarization_block10 = data10['start_depolarization_block']['evaluations']*1e-3

time_before_first_spike1 = data1['time_before_first_spike']['evaluations']
time_before_first_spike5 = data5['time_before_first_spike']['evaluations']
time_before_first_spike10 = data10['time_before_first_spike']['evaluations']

start_depolarization_block = [start_depolarization_block_nominal, start_depolarization_block1, start_depolarization_block5, start_depolarization_block10]
time_before_first_spike = [ time_before_first_spike_nominal, time_before_first_spike1, time_before_first_spike5, time_before_first_spike10]


# Standardized data
nominal = np.full(254,0)

start_depolarization_block1_n = zscore(start_depolarization_block1)
start_depolarization_block5_n = zscore(start_depolarization_block5)
start_depolarization_block10_n = zscore(start_depolarization_block10)

time_before_first_spike1_n = zscore(time_before_first_spike1)
time_before_first_spike5_n = zscore(time_before_first_spike5)
time_before_first_spike10_n = zscore(time_before_first_spike10)

start_depolarization_block_n = [nominal, start_depolarization_block1_n, start_depolarization_block5_n, start_depolarization_block10_n]
time_before_first_spike_n = [ nominal, time_before_first_spike1_n, time_before_first_spike5_n, time_before_first_spike10_n]

# Setup figure
set_style("seaborn-paper")
fig = plt.figure()
ax = fig.add_gridspec(3, 2)
axesA = fig.add_subplot(ax[0,0])
axesB = fig.add_subplot(ax[0,1])
axesC = fig.add_subplot(ax[1,0])
axesD = fig.add_subplot(ax[1,1])
axesE = fig.add_subplot(ax[2, 0:2])

fig.set_figwidth(9)
fig.set_figheight(9)

colors = ["#9467bd", "#17becf"]

# Panel A - B
axesA.hist(start_depolarization_block5, bins=15, color=colors[0])
axesA.axvline(start_depolarization_block_nominal[0], color='r', linestyle='dashed', linewidth=2)
axesB.hist(time_before_first_spike5, bins=15, color=colors[1])
axesB.axvline(time_before_first_spike_nominal[0], color='r', linestyle='dashed', linewidth=2)

axesA.set_xlabel(r"$T_\mathrm{sDP}$ $\mathrm{[s]}$", size = 13)
axesB.set_xlabel(r"$T_\mathrm{bFAP}$ $\mathrm{[ms]}$", size = 13)

axesA.set_yticks([])
axesB.set_yticks([])

axesA.text(3.43, 37, 'Baseline', fontsize=9)
axesB.text(3.68, 40, 'Baseline', fontsize=9)

# Violin plots
# Panel C
parts_nap = axesC.violinplot(
        start_depolarization_block, showmeans=False, showmedians=False,
        showextrema=False)
for pc in parts_nap['bodies']:
    pc.set_facecolor('#9467bd')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

quartile1, medians, quartile3 = np.percentile(start_depolarization_block, [25, 50, 75], axis=1)

inds = np.arange(1, len(medians) + 1)
axesC.scatter(inds[1:], medians[1:], marker='o', color='white', s=10, zorder=3)
axesC.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=3)

# Panel D
parts_tbfs = axesD.violinplot(
        time_before_first_spike, showmeans=False, showmedians=False,
        showextrema=False)
for pc in parts_tbfs['bodies']:
    pc.set_facecolor('#17becf')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

quartile1, medians, quartile3 = np.percentile(time_before_first_spike, [25, 50, 75], axis=1)

inds = np.arange(1, len(medians) + 1)
axesD.scatter(inds[1:], medians[1:], marker='o', color='white', s=10, zorder=3)
axesD.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=3)

# Panel E
# S_DP
quartile1, medians, quartile3 = np.percentile(start_depolarization_block_n, [25, 50, 75], axis=1)
inds = np.arange(1, len(medians) + 1)

parts_sdp_n = axesE.violinplot(
        start_depolarization_block_n,  positions=inds-0.11, widths=0.2, showmeans=False, showmedians=False,
        showextrema=False)

for pc in parts_sdp_n['bodies']:
    pc.set_facecolor('#9467bd')
    pc.set_edgecolor('black')
    pc.set_alpha(1)


axesE.scatter(inds[1:]-0.11, medians[1:], marker='o', color='white', s=10, zorder=3)
axesE.vlines(inds-0.11, quartile1, quartile3, color='k', linestyle='-', lw=3)

# T_bFS
quartile1, medians, quartile3 = np.percentile(time_before_first_spike_n, [25, 50, 75], axis=1)
inds = np.arange(1, len(medians) + 1)

parts_tbfs_n = axesE.violinplot(
        time_before_first_spike_n,  positions=inds+0.11, widths=0.2, showmeans=False, showmedians=False,
        showextrema=False)

for pc in parts_tbfs_n['bodies']:
    pc.set_facecolor('#17becf')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

axesE.scatter(inds[1:]+0.11, medians[1:], marker='o', color='white', s=10, zorder=3)
axesE.vlines(inds+0.11, quartile1, quartile3, color='k', linestyle='-', lw=3)

axesE.legend([parts_sdp_n['bodies'][0], parts_tbfs_n['bodies'][0]],[r"$T_\mathrm{sDP}^\mathrm{(std)}$", r"$T_\mathrm{bFAP}^\mathrm{(std)}$"], 
              fontsize=12 ,bbox_to_anchor=(0.2,1.2), frameon=False, labelspacing=0.1)

# set style for the axes
labels = ['0', '1', '5', '10']
for ax in [axesC, axesD, axesE]:
    set_axis_style(ax, labels)


axesC.set_ylabel(r"$T_\mathrm{sDP}$ $\mathrm{[s]}$", size=13)
axesD.set_ylabel(r"$T_\mathrm{bFAP}$ $\mathrm{[ms]}$", size=13)

# ABC
axesA.text(-0.05, 1.2, 'A', transform=axesA.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesB.text(-0.05, 1.2, 'B', transform=axesB.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesC.text(-0.05, 1.2, 'C', transform=axesC.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesD.text(-0.05, 1.2, 'D', transform=axesD.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')
axesE.text(-0.05, 1.2, 'E', transform=axesE.transAxes, fontsize=17, fontweight='bold', va='top', ha='right')

# Set subtitles
axesB.set_title(r"Histograms of QoIs ($\sigma$ = 5%)", fontsize=15, x=-0.2, y=1.3)
axesD.set_title(r"Violin plots of QoIs", fontsize=15, x=-0.2, y=1.3)
axesE.set_title(r"Grouped violin plot of standardized QoIs", fontsize=15, x=0.5, y=1.3)

plt.tight_layout()
plt.subplots_adjust(wspace=0.4, hspace=1, top=0.89)

# Save path
# checking if the directory exist and create it if it doesn't
if not os.path.exists('full_figures'):
    os.makedirs('full_figures')
plt.savefig('full_figures/Figure6.pdf', dpi=600)
# plt.show()