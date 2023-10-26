from __future__ import division, print_function
import bilby
from gwpy.timeseries import TimeSeries
import json
import os


outdir = 'outdir_GW190915_235702_GR'
label = 'GW190915_235702'
logger = bilby.core.utils.logger
bilby.core.utils.setup_logger(outdir=outdir, label=label)

datapath = os.path.expanduser('~/works/data/GW_data/data_info.json')
with open(datapath, 'r') as f:
    data_contents = json.load(f)
event_info = data_contents[label]

interferometer_names = list(event_info['strain'].keys())
time_of_event = event_info['gpstime']
start_time = time_of_event - 2
end_time = time_of_event + 2
duration = end_time - start_time
sampling_frequency = 2048
low_freq_dict = dict(L1=20., H1=20., V1=20.)        # low frequency cutoff referring the GWTC-1 paper

# READ IN THE DATA
ifo_list = bilby.gw.detector.InterferometerList([])
data_dict = event_info['strain']
psd_dict = event_info['PSDs']
for det in interferometer_names:
    logger.info("Reading analysis data for ifo {}".format(det))
    ifo = bilby.gw.detector.get_empty_interferometer(det)
    data = TimeSeries.read(data_dict[det], format='hdf5.gwosc', start=start_time, end=end_time)
    data = data.resample(rate=2048)
    ifo.set_strain_data_from_gwpy_timeseries(data)
    ifo.minimum_frequency = low_freq_dict[det]
    logger.info("Reading psd data for ifo {}".format(det))
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(psd_file=psd_dict[det])
    ifo_list.append(ifo)
logger.info("Saving data plots to {}".format(outdir))
bilby.core.utils.check_directory_exists_and_if_not_mkdir(outdir)
ifo_list.plot_data(outdir=outdir, label=label)

# CHOOSE PRIOR FILE
priors = bilby.gw.prior.BBHPriorDict(filename=os.path.expanduser('~/works/codes/ground_pca/AlignedSpin_high.prior'))
deltaT = 1.0
priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=time_of_event - deltaT / 2,
    maximum=time_of_event + deltaT / 2,
    name='geocent_time',
    latex_label='$t_c$',
    unit='$s$')
priors['chirp_mass'] = bilby.core.prior.Uniform(name='chirp_mass', minimum=10, maximum=50)
priors['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.1, maximum=1.0)
priors['luminosity_distance'] = bilby.gw.prior.UniformSourceFrame(name='luminosity_distance', minimum=100, maximum=5000)

# GENERATE WAVEFORM
waveform_arguments = dict(waveform_approximant='IMRPhenomXAS', 
                          reference_frequency=20., 
                          minimum_frequency=20.,
                          maximum_frequency=1024.,
                          catch_waveform_errors=True)
waveform_generator = bilby.gw.WaveformGenerator(duration=duration,
                                                sampling_frequency=sampling_frequency,
                                                frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
                                                parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
                                                waveform_arguments=waveform_arguments,)

# CHOOSE LIKELIHOOD FUNCTION
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers=ifo_list, 
    waveform_generator=waveform_generator,
    priors=priors,
    time_marginalization=False,
    distance_marginalization=False,
    phase_marginalization=True)

# RUN SAMPLER
result = bilby.run_sampler(
    likelihood=likelihood,
    priors=priors,
    outdir=outdir,
    label=label,
    sampler='pymultinest',
    npoints=2048,
    # walks=200,
    # nact=10,
    # maxmcmc=15000,
    # queue_size=72,
    dlogz=0.1,
    conversion_function=bilby.gw.conversion.generate_all_bbh_parameters)

# PLOT RESULT
plot_parameter_keys_corner = ['mass_1_source','mass_2_source','chirp_mass_source',
                              'chirp_mass', 'total_mass_source','mass_ratio',
                              'luminosity_distance', 'a_1','a_2','chi_eff']
result.plot_corner(parameters=plot_parameter_keys_corner, 
                   quantiles=[0.05, 0.95])


# COMPARISION
from pesummary.io import read
from pesummary.core.plots.plot import _make_comparison_corner_plot
import matplotlib
params = {'text.usetex': False, }
matplotlib.rcParams.update(params)

file_name_bilby = outdir+'/{}_result.json'.format(label)
data_bilby = read(file_name_bilby)
posterior_samples_bilby = data_bilby.samples_dict

file_name_LIGO = os.path.expanduser('~/works/data/GW_data/posterior_and_psd/GWTC-2_all_posterior_samples/{}.h5'.format(label))
data_LIGO = read(file_name_LIGO)
posterior_samples_LIGO = data_LIGO.samples_dict['PublicationSamples']
samples_dict = {'bilby':posterior_samples_bilby, 'LIGO_PosteriorSamples':posterior_samples_LIGO}

plot_parameter_keys_comparision = ['mass_1_source','mass_2_source','chirp_mass_source',
                                   'chirp_mass', 'total_mass_source','mass_ratio',
                                   'luminosity_distance', 'a_1','a_2','chi_eff']
fig = _make_comparison_corner_plot(samples_dict, latex_labels=posterior_samples_bilby.latex_labels, 
                                   corner_parameters=plot_parameter_keys_comparision, quantiles=[0.05, 0.95])
fig.savefig(outdir + '/comparision_corner_plot.png')
