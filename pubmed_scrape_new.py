import os
import time
import random
import datetime
import shutil

from Bio import Entrez

import util

# proj = 'terrawhisper'
# query = f'medicinal plants'.strip().lower()
proj = 'ozonogroup'
query = f'ozone'.strip().lower()
query_slug = query.replace(' ', '-')

vault = '/home/ubuntu/vault'
pubmed_folderpath = f'{vault}/{proj}/studies/pubmed'
query_folderpath = f'{pubmed_folderpath}/{query_slug}'

Entrez.email = 'martinpellizzer@gmail.com'
# sort_by = 'pub_date'
sort_by = 'relevance'
datetypes = ['mdat', 'pdat', 'edat']
datetype = datetypes[2]
years = [year for year in range(2025, 2020, -1)]
yesterday = datetime.datetime.now() - datetime.timedelta(1)
yesterday = datetime.datetime.strftime(yesterday, '%Y/%m/%d')

actions_num_total = 0

def create_folder(folderpath):
    chunk_curr = ''
    for chunk in folderpath.split('/'):
        chunk_curr += f'{chunk}/'
        try: os.makedirs(chunk_curr)
        except: continue

def get_ids(query, year=None, date=None, retmax=50):
    if date:
        handle = Entrez.esearch(db='pubmed', term=query, retmax=retmax, sort=sort_by, datetype=datetype, mindate=date, maxdate=date)
    elif year:
        handle = Entrez.esearch(db='pubmed', term=query, retmax=retmax, sort=sort_by, datetype=datetype, mindate=year, maxdate=year)
    else:
        handle = Entrez.esearch(db='pubmed', term=query, retmax=retmax, sort=sort_by)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def fetch_details(pmid):
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    try: records = Entrez.read(handle)
    except: 
        handle.close()
        return None
    handle.close()
    return records

def scrape_pubmed_jsons(query, year=None):
    global actions_num_total
    query_slug = query.strip().lower().replace(' ', '-')
    pmid_list = get_ids(query, year, retmax=9999)
    if pmid_list:
        create_folder(f'{query_folderpath}/json')
        for pmid_i, pmid in enumerate(pmid_list):
            done_pmids = []
            for filename in os.listdir(f'{query_folderpath}/json'):
                done_pmid = filename.split('.')[0]
                done_pmids.append(done_pmid)
            print(f'{pmid_i}/{len(pmid_list)} - {year}')
            if pmid in done_pmids: continue
            try: details = fetch_details(pmid)
            except: continue
            util.json_write(f'{query_folderpath}/json/{pmid}.json', details)
            time.sleep(random.randint(2, 5))
            actions_num_total += 1
            if actions_num_total >= 1000:
                actions_num_total = 0
                time.sleep(random.randint(450, 750))
    else:
        print('no article found')

for year in years:
    scrape_pubmed_jsons(query, year)
    
