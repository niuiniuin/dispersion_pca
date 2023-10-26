import numpy as np
import json
import os
import h5py


# GWTC1_events = ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 
#                 'GW170729', 'GW170809', 'GW170814', 'GW170817', 'GW170818', 
#                 'GW170823']
# GWTC2_events = ['GW190408_181802', 'GW190412',        'GW190413_052954', 'GW190413_134308', 
#                 'GW190421_213856', 'GW190424_180648', 'GW190425',        'GW190426_152155', 
#                 'GW190503_185404', 'GW190512_180714', 'GW190513_205428', 'GW190514_065416', 
#                 'GW190517_055101', 'GW190519_153544', 'GW190521',        'GW190521_074359', 
#                 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', 
#                 'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 
#                 'GW190719_215514', 'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 
#                 'GW190731_140936', 'GW190803_022701', 'GW190814',        'GW190828_063405', 
#                 'GW190828_065509', 'GW190909_114149', 'GW190910_112807', 'GW190915_235702', 
#                 'GW190924_021846', 'GW190929_012149', 'GW190930_133541']


# with open('data_info.json', 'r') as f:
#     data_info = json.load(f)

# with open('GWTC12.json', 'r') as f:
#     events_json = json.load(f)
#     events_json = events_json['events']
# GWTC1_events = []
# GWTC2_events = []
# for event, info in events_json.items():
#     catalog = info['catalog.shortName']
#     eventname = info['commonName']
#     if catalog=='GWTC-1-confident':
#         GWTC1_events.append(eventname)
#     elif catalog=='GWTC-2':
#         GWTC2_events.append(eventname)


# for event in GWTC1_events:
#     PSDs = np.loadtxt('GWTC1_PSDs/GWTC1_{}_PSDs.dat'.format(event))
#     if len(PSDs[0])==4:
#         H1_psd = PSDs[:,[0,1]]
#         L1_psd = PSDs[:,[0,2]]
#         V1_psd = PSDs[:,[0,3]]
#         np.savetxt('{}/{}_H1_psd.txt'.format(event, event), H1_psd)
#         np.savetxt('{}/{}_L1_psd.txt'.format(event, event), L1_psd)
#         np.savetxt('{}/{}_V1_psd.txt'.format(event, event), V1_psd)
#         data_info[event]['PSDs'] = {'H1':os.path.abspath('{}/{}_H1_psd.txt'.format(event, event)),
#                                     'L1':os.path.abspath('{}/{}_L1_psd.txt'.format(event, event)),
#                                     'V1':os.path.abspath('{}/{}_V1_psd.txt'.format(event, event))}
#     elif len(PSDs[0])==3:
#         H1_psd = PSDs[:,[0,1]]
#         L1_psd = PSDs[:,[0,2]]
#         np.savetxt('{}/{}_H1_psd.txt'.format(event, event), H1_psd)
#         np.savetxt('{}/{}_L1_psd.txt'.format(event, event), L1_psd)
#         data_info[event]['PSDs'] = {'H1':os.path.abspath('{}/{}_H1_psd.txt'.format(event, event)),
#                                     'L1':os.path.abspath('{}/{}_L1_psd.txt'.format(event, event))}

# for event in GWTC2_events:
#     filename = '/home/hydrogen/workspace/parity_new/LVC_posterior/all_posterior_samples/{}.h5'.format(event)
#     with h5py.File(filename, 'r') as f:
#         for k in f.keys():
#             if 'psds' in f[k].keys():
#                 PSDs = f[k]['psds']
#                 dets = list(PSDs.keys())
#                 if len(dets)>0:
#                     psds_dict = {}
#                     for detname in dets:
#                         det_psd = np.array(PSDs[detname])
#                         np.savetxt('{}/{}_{}_psd.txt'.format(event, event, detname), det_psd)
#                         psds_dict[detname] = os.path.abspath('{}/{}_{}_psd.txt'.format(event, event, detname))
#                     continue
#                 else:
#                     pass
#             else:
#                 pass
#     data_info[event]['PSDs'] = psds_dict


# with open('data_info.json', 'w') as f:
#     json.dump(data_info, f, indent=2)















# ['C01:IMRPhenomD', 'C01:IMRPhenomPv2', 'C01:NRSur7dq4', 'C01:SEOBNRv4P', 'C01:SEOBNRv4P_nonevol', 'PrecessingSpinIMR', 'PrecessingSpinIMRHM', 'PublicationSamples', 'ZeroSpinIMR', 'history', 'version']
filename = '/home/hydrogen/workspace/parity_new/LVC_posterior/all_posterior_samples/{}.h5'.format('GW190519_153544')
with h5py.File(filename, 'r') as f:
    print(f.keys())
    key = 'C01:IMRPhenomPv2'
    print(f[key]['psds'].keys())
    H1psd = np.array(f[key]['psds']['H1'])
    L1psd = np.array(f[key]['psds']['L1'])
    V1psd = np.array(f[key]['psds']['V1'])
print(H1psd[-10:])
print(L1psd[-10:])
print(H1psd[-10:])


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
