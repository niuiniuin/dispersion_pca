import bilby
from scipy.stats import gaussian_kde
import numpy as np
import os
import sys
sys.path.append(os.path.expanduser('~/works/codes/ground_pca'))
from  waveform_dispersion import waveform_alpha_all


'''
        commonName
0  GW190706_222641
1  GW190602_175927
2  GW190519_153544
3  GW190513_205428
4  GW190828_063405
5         GW170823
6  GW190915_235702
7  GW190408_181802
8  GW190503_185404
9  GW190512_180714
'''

dispersion_paras = ['A_alpha_00', 'A_alpha_25', 'A_alpha_30', 'A_alpha_40']
posterior_1  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190513_205428_dispersion/GW190513_205428_result.json').posterior[dispersion_paras]
posterior_2  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190519_153544_dispersion/GW190519_153544_result.json').posterior[dispersion_paras]
posterior_3  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190602_175927_dispersion/GW190602_175927_result.json').posterior[dispersion_paras]
posterior_4  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190706_222641_dispersion/GW190706_222641_result.json').posterior[dispersion_paras]
posterior_5  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190828_063405_dispersion/GW190828_063405_result.json').posterior[dispersion_paras]
posterior_6  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW170823_dispersion/GW170823_result.json').posterior[dispersion_paras]
posterior_7  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190915_235702_dispersion/GW190915_235702_result.json').posterior[dispersion_paras]
posterior_8  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190408_181802_dispersion/GW190408_181802_result.json').posterior[dispersion_paras]
posterior_9  = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190503_185404_dispersion/GW190503_185404_result.json').posterior[dispersion_paras]
posterior_10 = bilby.read_in_result('/home/changfenggroup/nrui/works/codes/ground_pca/runs/outdir_GW190512_180714_dispersion/GW190512_180714_result.json').posterior[dispersion_paras]

# print(posterior_1.transpose().shape)
kde_1  = gaussian_kde(posterior_1.transpose())
kde_2  = gaussian_kde(posterior_2.transpose())
kde_3  = gaussian_kde(posterior_3.transpose())
kde_4  = gaussian_kde(posterior_4.transpose())
kde_5  = gaussian_kde(posterior_5.transpose())
kde_6  = gaussian_kde(posterior_6.transpose())
kde_7  = gaussian_kde(posterior_7.transpose())
kde_8  = gaussian_kde(posterior_8.transpose())
kde_9  = gaussian_kde(posterior_9.transpose())
kde_10 = gaussian_kde(posterior_10.transpose())

class Likelihood_KDE(bilby.core.likelihood.Likelihood):

    def __init__(self):
        super(Likelihood_KDE, self).__init__(parameters={})

    def log_likelihood(self):
        return (
                kde_1.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_2.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_3.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_4.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_5.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_6.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_7.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_8.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
                kde_9.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']]) +
               kde_10.logpdf([self.parameters['A_alpha_00'], self.parameters['A_alpha_25'], self.parameters['A_alpha_30'], self.parameters['A_alpha_40']])
                )

priors = {}
priors['A_alpha_00'] = bilby.core.prior.Uniform(name='A_alpha_00', minimum=-5.0, maximum=5.0)
priors['A_alpha_25'] = bilby.core.prior.Uniform(name='A_alpha_25', minimum=-300.0, maximum=300.0)
priors['A_alpha_30'] = bilby.core.prior.Uniform(name='A_alpha_30', minimum=-300.0, maximum=300.0)
priors['A_alpha_40'] = bilby.core.prior.Uniform(name='A_alpha_40', minimum=-50.0, maximum=50.0)

likelihood = Likelihood_KDE()
# likelihood.parameters.update({'A_alpha_00':0.99931872, 'A_alpha_25':-13.62101078, 'A_alpha_30':-17.08741665, 'A_alpha_40':-3.58600163})
# print(likelihood.log_likelihood())

# RUN SAMPLER
label = 'combining'
outdir = 'outdir_'+label
result = bilby.run_sampler(
    likelihood=likelihood,
    priors=priors,
    outdir=outdir,
    label=label,
    nlives=1024,
    sampler='pymultinest')

# PLOT RESULT
result.plot_corner(quantiles=[0.05, 0.95],
                   parameters={'A_alpha_00': 0.0, 
                               'A_alpha_25': 0.0, 
                               'A_alpha_30': 0.0, 
                               'A_alpha_40': 0.0}
                   )
