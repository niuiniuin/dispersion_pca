import matplotlib 
from matplotlib import pyplot as plt
import bilby
from numpy.linalg import eig, inv
import numpy as np
import pandas as pd
import scipy.linalg as sci_la

fig_width_pt = 3*246.0                  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0      # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = { 'axes.labelsize': 35,
          'font.family': 'serif',
          'font.serif': 'Computer Modern Raman',
          'font.size': 35,
          'legend.fontsize': 20,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'axes.grid' : True,
          'text.usetex': True,
          'savefig.dpi' : 100,
          'lines.markersize' : 6,
          'figure.figsize': fig_size}
matplotlib.rcParams.update(params)

########################################################################################################################
# import numpy as np
# from numpy.linalg import eig, inv
# A = np.array([[1,4,8],[3,5,3],[3,3,1]])
# A_eval, A_evec = eig(A)
# # print(A_eval)
# # print(A_evec)
# idx = A_eval.argsort()
# A_eval = A_eval[idx]
# A_evec = A_evec[:,idx]
# # print(A_eval)
# # print(A_evec)

# # print(A @ A_evec[:,0].reshape((3,1)))
# # print(A_eval[0] * A_evec[:,0].reshape((3,1)))
# # print(A @ A_evec[:,1].reshape((3,1)))
# # print(A_eval[1] * A_evec[:,1].reshape((3,1)))
# # print(A @ A_evec[:,2].reshape((3,1)))
# # print(A_eval[2] * A_evec[:,2].reshape((3,1)))

# # print(A_eval)
# # print(A_evec)
# A_eval_mat = np.diag(A_eval)
# # print(A_eval_mat)
# A_evec_inv = inv(A_evec)

# print(A - (A_evec@A_eval_mat@A_evec_inv))
# print((A@A_evec) - (A_evec@A_eval_mat))

# print(A_evec@A_eval_mat@A_evec_inv)
# print(A@A_evec)
# print(A_evec@A_eval_mat)
########################################################################################################################
res_dir = 'outdir_combining'
res_filename = res_dir + '/combining_result.json'
res_file = bilby.result.read_in_result(filename=res_filename)
# print(res_file.label)
res_label = res_file.label
res_posterior = res_file.posterior
# print(res_posterior)
tgr_para = ['A_alpha_00', 'A_alpha_25', 'A_alpha_30', 'A_alpha_40']
tgr_para_posterior = res_posterior[tgr_para]

# compute covariance matrix, and eigen vector eigen value of it
cov_tgr_para = tgr_para_posterior.cov()
# print(cov_tgr_para)
tgr_eval, tgr_evec = eig(cov_tgr_para)
idx = tgr_eval.argsort()
tgr_eval_sorted = tgr_eval[idx]
tgr_evec_sorted = tgr_evec[:,idx]
# print('eval:', tgr_eval_sorted)
# # print for check
# print('cov matrix of tgr_para: \n', cov_tgr_para)
# print('eigen values: \n', tgr_eval)
# print('eigen evectors: \n', tgr_evec)
# print('sorted eigen values: \n', tgr_eval_sorted)
# print('sorted eigen vectors: \n', tgr_evec_sorted)
# print(tgr_evec_sorted[:,0])
# print(tgr_evec_sorted[:,0]@tgr_evec_sorted[:,0])

# tgr_eval_sci, tgr_evec_sci = eig(cov_tgr_para)
# idx_sci = tgr_eval_sci.argsort()
# tgr_eval_sorted_sci = tgr_eval_sci[idx_sci]
# tgr_evec_sorted_sci = tgr_evec_sci[:,idx_sci]
# print(tgr_eval_sorted_sci)
# print(tgr_evec_sorted_sci[:,0].tgr_evec_sorted_sci[:,0])

# # for test
# print('inv, T', np.sum((tgr_evec_sorted.T-inv(tgr_evec_sorted)).flatten()))
# print('A v=lambda v', np.sum((np.array(cov_tgr_para)@tgr_evec - tgr_evec@np.diag(tgr_eval)).flatten()))
# print('A = vT lambda v', np.sum((np.array(cov_tgr_para)- (tgr_evec@np.diag(tgr_eval)@tgr_evec.T)).flatten()))
# print('A = v-1 lambda v', np.sum((np.array(cov_tgr_para)- (tgr_evec@np.diag(tgr_eval)@inv(tgr_evec))).flatten()))
########################################################################################################################
tgr_posteror_PCA = tgr_para_posterior@tgr_evec_sorted
# print(tgr_para_posterior)
# print(tgr_posteror_PCA)
# print(tgr_posteror_PCA[0])
PCA_para_name = [f'PCA_{i}' for i in range(len(tgr_para))]
for i in range(len(tgr_para)):
    print(f'{PCA_para_name[i]} and {tgr_para[i]}')
    res_posterior['PCA_{}'.format(i)] = tgr_posteror_PCA[i]
    # print(f'standard deviation of {i}th principal component: ', np.std(tgr_posteror_PCA[i]))
    # print(f'standard deviation of {i}th tgr_para: ', np.std(tgr_para_posterior[tgr_para[i]]))

    lower, median, upper = np.percentile(tgr_posteror_PCA[i], [0.15, 50, 99.85])
    sigma = np.std(tgr_posteror_PCA[i])
    print(f'{i} pca: 3-sigma bound: ', lower, upper)
    print(f'{i} pca: GR value at ', (0.0 - median)/sigma)
    lower, median, upper = np.percentile(tgr_para_posterior[tgr_para[i]], [0.15, 50, 99.85])
    sigma = np.std(tgr_para_posterior[tgr_para[i]])
    print(f'{i} original: 3-sigma bound: ', lower, upper)
    print(f'{i} original: GR value at ', (0.0 - median)/sigma)

    # print(f'percentile 0.05 and 0.95 of {i}th principal component: ', np.percentile(tgr_posteror_PCA[i], [5, 95]))
    # print(f'percentile 0.05 and 0.95 of {i}th tgr_para: ', np.percentile(tgr_para_posterior[tgr_para[i]], [5, 95]))

# dict_tgr_para = dict.fromkeys(tgr_para, 0.0)
# dict_PCA_para_name = dict.fromkeys(PCA_para_name, 0.0)
# print(dict_tgr_para, dict_PCA_para_name)
dict_dispersion_para = {}
# dict_dispersion_para.update({"A_alpha_00": 0.005,
#                              "A_alpha_05": 0.03,
#                              "A_alpha_10": 0.2,
#                              "A_alpha_15": 0.3,
#                              "A_alpha_25": 1.2,
#                              "A_alpha_30": 1.8,
#                              "A_alpha_35": 1.6,
#                              "A_alpha_40": 0.2,})
dict_dispersion_para.update({"A_alpha_00": 0.0,
                             "A_alpha_25": 0.0,
                             "A_alpha_30": 0.0,
                             "A_alpha_40": 0.0,})
inj_para_pd = pd.DataFrame([dict_dispersion_para])
PCA_inj_para = inj_para_pd@tgr_evec_sorted
dict_PCA_para = {}
for idx, pca_para in enumerate(PCA_para_name):
    dict_PCA_para[pca_para] = PCA_inj_para.loc[0, idx]
# print(dict_tgr_para)
# print(dict_PCA_para)


tgr_para_labels = [r'$\delta\phi_{0.0}$',
                   r'$\delta\phi_{2.5}$', 
                   r'$\delta\phi_{3.0}$', 
                   r'$\delta\phi_{4.0}$']
PCA_para_labels = [r'$\delta\phi^{\rm PCA}_{0\,\rm th}$', 
                   r'$\delta\phi^{\rm PCA}_{1\,\rm st}$', 
                   r'$\delta\phi^{\rm PCA}_{2\,\rm nd}$', 
                   r'$\delta\phi^{\rm PCA}_{3\,\rm rd}$']
res_file.plot_corner(parameters=dict_dispersion_para, quantiles=[0.0015, 0.9985], titles=True, labels=tgr_para_labels, label_kwargs=dict(fontsize=24), filename=f'{res_dir}/dispersion_{res_label}.pdf')
res_file.plot_corner(parameters=dict_PCA_para, quantiles=[0.0015, 0.9985], titles=True, labels=PCA_para_labels, label_kwargs=dict(fontsize=24), filename=f'{res_dir}/dispersion_{res_label}_pca.pdf')
