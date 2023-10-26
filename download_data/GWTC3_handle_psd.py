import numpy as np
import json
import os
import h5py



with open('../data_info.json', 'r') as f:
    data_info = json.load(f)

with open('GWTC3.json', 'r') as f:
    events_json = json.load(f)
    events_json = events_json['events']
GWTC3_events = []
for event, info in events_json.items():
    catalog = info['catalog.shortName']
    eventname = info['commonName']
    if catalog=='GWTC-3-confident':
        GWTC3_events.append(eventname)
print(GWTC3_events)
print(len(GWTC3_events))



for event in GWTC3_events:
    psds_dict = {}
    filename = '/home/hydrogen/workspace/GW_data/posterior_and_psd/GWTC-3_posterior/IGWN-GWTC3p0-v1-{}_PEDataRelease_mixed_nocosmo.h5'.format(event)
    with h5py.File(filename, 'r') as f:
        for k in f.keys():
            if 'psds' in f[k].keys():
                PSDs = f[k]['psds']
                dets = list(PSDs.keys())
                if len(dets)>0:
                    print('find the psd data of event: {}, stored in the label: {}, contains dets: {}'.format(event, k, dets))
                    for detname in dets:
                        det_psd = np.array(PSDs[detname])
                        np.savetxt('/home/hydrogen/workspace/GW_data/{}/{}_{}_psd.txt'.format(event, event, detname), det_psd)
                        psds_dict[detname] = '/home/hydrogen/workspace/GW_data/{}/{}_{}_psd.txt'.format(event, event, detname)
                    # continue next event
                    break
        if psds_dict==0:
            print('cannot find psd data for event: {}'.format(event))

    data_info[event]['PSDs'] = psds_dict


with open('../data_info.json', 'w') as f:
    json.dump(data_info, f, indent=2)














# filename = '/home/hydrogen/workspace/GW_data/posterior_and_psd/GWTC-3_posterior/IGWN-GWTC3p0-v1-GW191103_012549_PEDataRelease_mixed_cosmo.h5'
# with h5py.File(filename, 'r') as f:
#     print(f.keys())
#     key = 'C01:IMRPhenomXPHM'
#     print(f[key]['psds'].keys())
#     H1psd = np.array(f[key]['psds']['H1'])
#     L1psd = np.array(f[key]['psds']['L1'])
#     V1psd = np.array(f[key]['psds']['V1'])
# print(H1psd[-10:])
# print(L1psd[-10:])
# print(H1psd[-10:])


    # for k in f.keys():
    #    if 'psds' in f[k].keys():
    #         PSDs = f[k]['psds']
    #         dets = list(PSDs.keys())
    #         if len(dets)>0:
    #                 psds_dict = {}
    #                 for detname in dets:
    #                     det_psd = np.array(PSDs[detname])
    #                     np.savetxt('{}/{}_{}_psd.txt'.format(event, event, detname), det_psd)
    #                     psds_dict[detname] = os.path.abspath('{}/{}_{}_psd.txt'.format(event, event, detname))




# ########################################################################################################################
# # the difference between cosmo and nocosmo
# from pesummary.io import read
# from pesummary.core.plots.plot import _make_comparison_corner_plot

# filename_cosmo   = '/home/hydrogen/workspace/GW_data/posterior_and_psd/GWTC-3_posterior/IGWN-GWTC3p0-v1-GW191103_012549_PEDataRelease_mixed_cosmo.h5'
# filename_nocosmo = '/home/hydrogen/workspace/GW_data/posterior_and_psd/GWTC-3_posterior/IGWN-GWTC3p0-v1-GW191103_012549_PEDataRelease_mixed_nocosmo.h5'

# data_cosmo = read(filename_cosmo)
# posterior_samples_cosmo = data_cosmo.samples_dict['C01:Mixed']

# data_nocosmo = read(filename_nocosmo)
# posterior_samples_nocosmo = data_nocosmo.samples_dict['C01:Mixed']

# plot_parameter_keys = ['mass_1_source','mass_2_source','chirp_mass_source',
#                        'chirp_mass', 'total_mass_source','mass_ratio',
#                        'luminosity_distance', 'a_1','a_2','chi_eff']
# comparision_dict = {'cosmo':posterior_samples_cosmo, 'nocosmo':posterior_samples_nocosmo}
# fig_comparision = _make_comparison_corner_plot(comparision_dict, latex_labels=posterior_samples_cosmo.latex_labels, 
#                                                corner_parameters=plot_parameter_keys, quantiles=[0.05, 0.95])
# fig_comparision.savefig('cosmo_nocosmo.png')
# fig_comparision.close()
# ########################################################################################################################
# for ename in GWTC3_events:
#     filename = '/home/hydrogen/workspace/GW_data/posterior_and_psd/GWTC-3_posterior/IGWN-GWTC3p0-v1-{}_PEDataRelease_mixed_cosmo.h5'.format(ename)
#     with h5py.File(filename, 'r') as f:
#         print(f.keys())
# ['C01:IMRPhenomXPHM', 'C01:Mixed', 'C01:SEOBNRv4PHM', 'history', 'version']
# ['C01:IMRPhenomNSBH:HighSpin', 'C01:IMRPhenomNSBH:LowSpin', 'C01:IMRPhenomXPHM:HighSpin', 'C01:IMRPhenomXPHM:LowSpin', 'C01:Mixed', 'C01:Mixed:NSBH:HighSpin', 'C01:Mixed:NSBH:LowSpin', 'C01:SEOBNRv4PHM', 'C01:SEOBNRv4_ROM_NRTidalv2_NSBH:HighSpin', 'C01:SEOBNRv4_ROM_NRTidalv2_NSBH:LowSpin', 'history', 'version']