import requests
import pandas as pd

# # download information of GWTC3
# gwtc_url = 'https://gwosc.org/eventapi/json/GWTC/'
# gwtc_json_source = urllib.request.urlopen(gwtc_url)
# print(gwtc_json_source.__dir__())
# with open('gwtc.json', 'bw') as f:
#     f.write(gwtc_json_source.read())


# select events for analysis
# criteria:
# 1. FAR < 10^3 y^-1
# 2. SNR > 12
# 3. p_astro > 0.99
# 4. the most far awary 5 events
url_requests = requests.get('https://gwosc.org/eventapi/json/GWTC/')
gwtc_info = (url_requests.json())['events']
keys_name = ['commonName', 'luminosity_distance', 'network_matched_filter_snr', 'far', 'p_astro']
selected_events = []
for name, data in gwtc_info.items():
    far = data['far']
    snr = data['network_matched_filter_snr']
    p_astro = data['p_astro']
    if (far < 1e-3) and (snr > 12) and (p_astro > 0.99):
        selected_events.append(data)

selected_events = pd.DataFrame(selected_events)[keys_name]
selected_events.sort_values('luminosity_distance', ascending=False, inplace=True, ignore_index=True)
# print(selected_events.iloc[:20])
print(selected_events)

'''
         commonName  luminosity_distance  network_matched_filter_snr           far  p_astro
0   GW190706_222641               3630.0                        13.4  5.000000e-05  1.00000
1   GW190602_175927               2840.0                        13.2  1.000000e-05  1.00000
2   GW190519_153544               2600.0                        15.9  1.000000e-05  1.00000
3   GW190513_205428               2210.0                        12.5  1.300000e-05  0.99998
4   GW190828_063405               2070.0                        16.5  1.000000e-05  1.00000
5          GW170823               1940.0                        12.2  1.000000e-07  1.00000
6   GW190915_235702               1750.0                        13.1  1.000000e-05  1.00000
7   GW190408_181802               1540.0                        14.6  1.000000e-05  1.00000
8   GW190503_185404               1520.0                        12.2  1.000000e-05  1.00000
9   GW190512_180714               1460.0                        12.7  1.000000e-05  1.00000

10  GW190521_074359               1080.0                        25.9  1.000000e-05  1.00000
11         GW170809               1030.0                        12.8  1.000000e-07  1.00000
12         GW170104                990.0                        13.8  1.000000e-07  1.00000
13  GW190708_232457                930.0                        13.4  3.100000e-04  0.99965
14  GW190728_064510                880.0                        13.1  1.000000e-05  1.00000
15  GW190630_185205                870.0                        16.4  1.000000e-05  1.00000
16  GW190707_093326                850.0                        13.1  1.000000e-05  1.00000
17         GW190412                720.0                        19.8  1.000000e-05  1.00000
18         GW170814                600.0                        17.7  1.000000e-07  1.00000
19         GW151226                450.0                        13.1  1.000000e-07  1.00000

20         GW150914                440.0                        26.0  1.000000e-07  1.00000
21         GW170608                320.0                        15.4  1.000000e-07  1.00000
22         GW190814                230.0                        25.3  1.000000e-05  1.00000
23         GW170817                 40.0                        33.0  1.000000e-07  1.00000
'''

