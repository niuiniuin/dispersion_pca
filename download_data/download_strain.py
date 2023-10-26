import json
import urllib.request
import os 
import sys
import math



'''handle the events with glitch substraction by hand'''
# event_glitch = ['GW170817', 
#                 'GW190413_134308', 'GW190424_180648', 'GW190503_185404', 'GW190513_205428', 
#                 'GW190514_065416', 'GW190701_203306', 'GW190924_021846',
#                 ]

def _progress(block_num, block_size, total_size):
    total_block = math.ceil(total_size/block_size)
    sys.stdout.write('\r{}/{}'.format(block_num, total_block))
    sys.stdout.flush()
    return None

# read in stored data information
with open('data_info.json', 'r') as f:
    data_info = json.load(f)

# read in the json file for GWTC events
with open('GWTC3.json', 'r') as f:
    events_json = json.load(f)
    events_json = events_json['events']

# download data for every event and generate the contents
events_to_download_list = []
events_existed_list = []

for event, info in events_json.items():
    ename = info['commonName']
    if ename in list(data_info.keys()):
        '''data of event has already exist'''
        events_existed_list.append(event)
    else:
        '''data of event to download'''
        events_to_download_list.append(event)

print('total num of events in catalog: {}'.format(len(list(events_json.keys()))))
print('num of events has already exist: {}'.format(len(events_existed_list)))
print('num of events to be download: {}'.format(len(events_to_download_list)))

events_dict = {}
for event in events_to_download_list:
    with urllib.request.urlopen(events_json[event]['jsonurl']) as response:
        event_details = json.load(response)
    event_details = event_details['events'][event]
    eventtime = event_details['GPS']
    eventname = event_details['commonName']
    straininfo = event_details['strain']
    straininfo = list(filter(lambda x: x['duration']==4096 and x['format']=='hdf5' and x['sampling_rate']==4096, straininfo))
    if os.path.exists(eventname):
        print('folder of {} are already exist, please check the code'.format(eventname))
        continue
    else:
        os.mkdir(eventname)
        strain_data = {}
        print('downloading strain data for event {}'.format(eventname))
        for item in straininfo:
            det = item['detector']
            print(det+': '+item['url'].split('/')[-1])
            filename, headers = urllib.request.urlretrieve(item['url'], eventname+'/'+item['url'].split('/')[-1], _progress)
            sys.stdout.write('\n')
            fullpath = os.path.abspath(filename)
            strain_data[det] = fullpath
        eventdict = {'gpstime': eventtime,
                     'strain': strain_data}
        events_dict[eventname] = eventdict

# store the contents
data_info.update(events_dict)
with open('data_info.json', 'w') as f:
    json.dump(data_info, f, indent=2)








# filename, headers = urllib.request.urlretrieve('https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5', 'H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5', _progress)
# print(filename)
# print(headers)
# print(os.path.abspath(filename))







# with urllib.request.urlopen('https://www.gw-openscience.org/eventapi/json/GWTC-1-confident/GW150914/v3') as res:
#     print(json.load(res)['events'].keys())