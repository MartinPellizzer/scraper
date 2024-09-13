import os
import time
import metapub
from urllib.request import urlretrieve
import textract
import json

from oliark import file_write

proj = 'ozonogroup'
vault = '/home/ubuntu/vault'
proj_path = f'{vault}/{proj}'
studies_path = f'{proj_path}/studies'
pubmed_path = f'{studies_path}/pubmed'
query_path = f'{pubmed_path}/ozone'
json_path = f'{query_path}/json'
pdfs_texts_done_path = f'{query_path}/pdfs-texts-done'
pdfs_texts_fail_path = f'{query_path}/pdfs-texts-fail'

pmids = [filename.split('.')[0] for filename in os.listdir(json_path)]
pmids_done = [filename.split('.')[0] for filename in os.listdir(pdfs_texts_done_path)]
pmids_fail = [filename.split('.')[0] for filename in os.listdir(pdfs_texts_fail_path)]

tot_num = 0
for i, pmid in enumerate(pmids):
    print(f'{i+1}/{len(pmids)} - {pmid}')
    if pmid in pmids_done: continue
    if pmid in pmids_fail: continue
    try: 
        src = metapub.FindIt(pmid)
    except:
        file_write(f'{pdfs_texts_fail_path}/{pmid}.txt', '')
        continue
    print(src.pma.title)
    # print(src.pma.abstract)
    if src.url:
        try:
            urlretrieve(src.url, 'test-path/f')
            text = textract.process(
                'test-path/f',
                extension='pdf',
                methodd='pdftotext',
                encoding='utf_8',
            )
            text = text.decode('utf-8')
            file_write(f'{pdfs_texts_done_path}/{pmid}.txt', text)
        except:
            file_write(f'{pdfs_texts_fail_path}/{pmid}.txt', '')
    else:
        file_write(f'{pdfs_texts_fail_path}/{pmid}.txt', '')
        print(src.reason)
    print()

    time.sleep(1)
    tot_num += 1
    if tot_num >= 1000:
        tot_num = 0
        time.sleep(600)

