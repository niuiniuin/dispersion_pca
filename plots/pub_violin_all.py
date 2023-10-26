
import matplotlib 
from matplotlib import pyplot as plt
import bilby
from numpy.linalg import eig, inv
import numpy as np
import scipy.linalg as sci_la
import pandas as pd
import os
import sys
sys.path.append(os.path.expanduser('~/works/codes/ground_pca'))
from  waveform_dispersion import waveform_alpha_all

fig_width_pt = 3*246.0                  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width*2,fig_height]
params = { 'axes.labelsize': 30,
          'font.family': 'serif',
          'font.serif': 'Computer Modern Raman',
          'font.size': 30,
          'legend.fontsize': 20,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'axes.grid' : False,
          'text.usetex': True,
          'savefig.dpi' : 100,
          'lines.markersize' : 14,
          'figure.figsize': fig_size}
matplotlib.rcParams.update(params)

########################################################################################################################
para_labels = [r'$\delta\phi_{0.0}$',
               r'$\delta\phi_{2.5}$', 
               r'$\delta\phi_{3.0}$', 
               r'$\delta\phi_{4.0}$']
GR_values = {"A_alpha_00": 0.0,
             "A_alpha_25": 0.0,
             "A_alpha_30": 0.0,
             "A_alpha_40": 0.0}
dispersion_paras = ['A_alpha_00', 'A_alpha_25', 'A_alpha_30', 'A_alpha_40']
events = ['GW190706_222641', 'GW190602_175927', 'GW190519_153544', 'GW190513_205428', 'GW190828_063405', 'GW170823', 'GW190915_235702', 'GW190408_181802', 'GW190503_185404', 'GW190512_180714']

for ename in events:
    res_filename = f'/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_{ename}_dispersion/{ename}_result.json'
    res_file = bilby.result.read_in_result(filename=res_filename)
    res_posterior = res_file.posterior[dispersion_paras]
    
    fig, ax_dis = plt.subplots(nrows=1, ncols=4, figsize=[fig_width, fig_height*0.8],subplot_kw={'box_aspect':3.0})
    
    for idx, para in enumerate(dispersion_paras):
        ax_dis[idx].violinplot(res_posterior[para], widths=1.2, showextrema=False, showmedians=False, quantiles={0.05, 0.95})
        lower, median, upper = np.percentile(res_posterior[para], [5, 50, 95])
        ylim = 1.1 * np.max(np.abs([lower, upper]))
        sigma = np.std(res_posterior[para])
        # print(sigma)
        # print(f'{idx} {para} GR value at ', (0.0 - median)/sigma)
        ax_dis[idx].vlines(1, lower, upper)
        ax_dis[idx].axhline(0.0, color='tab:red', linestyle='dashed', linewidth=1.8)
        ax_dis[idx].set_xticks([])
        ax_dis[idx].set_ylim(-ylim, ylim)
        ax_dis[idx].set_xlabel(para_labels[idx])

    fig.suptitle(ename, y=0.9)

    fig.tight_layout()
    fig.savefig(f'violin_{ename}.pdf')
    fig.clear()


