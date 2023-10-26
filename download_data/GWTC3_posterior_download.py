import json
import urllib.request
import os
import sys
import math
import hashlib


def _progress(block_num, block_size, total_size):
    total_block = math.ceil(total_size/block_size)
    sys.stdout.write('\r{}/{}'.format(block_num, total_block))
    sys.stdout.flush()
    return None



with open('GWTC3_posterior_download.json', 'r') as f:
    json_file = json.load(f)
    all_files = json_file['files']


for f in all_files:
    file_url = f['links']['self']
    file_name = f['key']
    file_name = os.path.join('GWTC-3_posterior', file_name)
    checksum = f['checksum'].split(':')[1]
    
    if os.path.exists(file_name):
        with open(file_name, 'rb') as existed_f:
            contents = existed_f.read()
            md5 = hashlib.md5(contents).hexdigest()
        if md5 == checksum:
            print(file_name, 'has already been downloaded')
            continue
    
    print('downloading {}'.format(file_name))
    filename, headers = urllib.request.urlretrieve(file_url, file_name, _progress)
    sys.stdout.write('\n')

print('download sucessful')