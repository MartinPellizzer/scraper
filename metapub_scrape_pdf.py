import os
import time
import metapub
from urllib.request import urlretrieve
import json

from oliark import file_write

# proj = 'ozonogroup'
# query = 'ozone'
proj = 'terrawhisper'
query = 'medicinal-plants'

vault = '/home/ubuntu/vault'
proj_path = f'{vault}/{proj}'
studies_path = f'{proj_path}/studies'
pubmed_path = f'{studies_path}/pubmed'
query_path = f'{pubmed_path}/{query}'
json_path = f'{query_path}/json'
pdfs_done_path = f'{query_path}/pdfs-done'
pdfs_fail_path = f'{query_path}/pdfs-fail'

try: os.makedirs(f'{pdfs_done_path}')
except: pass
try: os.makedirs(f'{pdfs_fail_path}')
except: pass

pmids = [filename.split('.')[0] for filename in os.listdir(json_path)]
pmids_done = [filename.split('.')[0] for filename in os.listdir(pdfs_done_path)]
pmids_fail = [filename.split('.')[0] for filename in os.listdir(pdfs_fail_path)]

tot_num = 0
for i, pmid in enumerate(pmids):
    print(f'{i+1}/{len(pmids)} - {pmid}')
    if pmid in pmids_done: continue
    if pmid in pmids_fail: continue
    try: 
        src = metapub.FindIt(pmid)
    except:
        file_write(f'{pdfs_fail_path}/{pmid}', '')
        continue
    print(src.pma.title)
    if src.url:
        try:
            urlretrieve(src.url, f'{pdfs_done_path}/{pmid}.pdf')
        except:
            file_write(f'{pdfs_fail_path}/{pmid}.pdf', '')
    else:
        file_write(f'{pdfs_fail_path}/{pmid}.txt', '')
        print(src.reason)
    print()

    time.sleep(1)
    tot_num += 1
    if tot_num >= 1000:
        tot_num = 0
        time.sleep(600)

