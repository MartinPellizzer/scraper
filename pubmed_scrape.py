import os
import time
import random
import datetime
import shutil

from Bio import Entrez

import util
import util_ai

model = '/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf'

# proj = 'terrawhisper'
proj = 'ozonogroup'
# query = f'medicinal plants'.strip().lower()
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
years = [year for year in range(2025, 1810, -1)]
yesterday = datetime.datetime.now() - datetime.timedelta(1)
yesterday = datetime.datetime.strftime(yesterday, '%Y/%m/%d')

actions_num_total = 0

def create_folder(folderpath):
    chunk_curr = ''
    for chunk in folderpath.split('/'):
        chunk_curr += f'{chunk}/'
        try: os.makedirs(chunk_curr)
        except: continue

def get_ids_yesterday(query, date=None, retmax=50):
    handle = Entrez.esearch(db='pubmed', term=query, retmax=retmax, sort=sort_by, datetype=datetype, mindate=date, maxdate=date)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def get_ids_latest(query, retmax=50):
    handle = Entrez.esearch(db='pubmed', term=query, retmax=retmax, sort='pub_date')
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

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

def scrape_pubmed_yesterday(query, yesterday):
    query_slug = query.strip().lower().replace(' ', '-')
    pmid_list = get_ids_yesterday(query, date=yesterday, retmax=50)
    # pmid_list = get_ids_latest(query, retmax=50)
    citation_arr = []
    abstract_arr = []
    try: shutil.rmtree(f'{query_folderpath}/yesterday')
    except: pass
    try: os.makedirs(f'{query_folderpath}/yesterday')
    except: pass
    shutil.rmtree(f'{query_folderpath}/yesterday.csv')
    row = [
        'pmid',
        'pub_year',
        'pub_month',
        'pub_day',
        'journal_title',
        'journal_volume',
        'journal_issue',
        'pages',
        'title',
        'author_st',
    ]
    util.csv_add_rows(f'{query_folderpath}/yesterday.csv', [row])
    if pmid_list:
        create_folder(query_folderpath)
        for pmid_i, pmid in enumerate(pmid_list):
            print(f'{pmid_i}/{len(pmid_list)}')
            abstract_text = ''
            title = ''
            authors = ''
            authors_str = ''
            journal_title = ''
            journal_volume = ''
            journal_issue = ''
            pub_date = ''
            pub_year = ''
            pub_month = ''
            pub_day = ''
            pages = ''
            try: 
                details = fetch_details(pmid)
            except:
                continue
            try:
                abstract_text = details['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            except:
                continue
            article = details['PubmedArticle'][0]['MedlineCitation']['Article']
            journal = article['Journal']
            pubmed_data = details['PubmedArticle'][0]['PubmedData']
            title = article.get('ArticleTitle', 'No title available')
            try: authors = article['AuthorList']
            except: pass
            try: author_str = ', '.join([f"{a['LastName']} {a['ForeName'][0]}" for a in authors])
            except: pass
            try: journal_title = journal.get('Title', 'No journal title available')
            except: pass
            try: journal_volume = journal['JournalIssue'].get('Volume', 'No volume')
            except: pass
            try: journal_issue = journal['JournalIssue'].get('Issue', 'No issue')
            except: pass
            try: pub_date = article.get('ArticleDate', [{'Year': 'No year', 'Month': 'No month', 'Day': 'No day'}])[0]
            except: pass
            try: pub_year = pub_date['Year']
            except: pass
            try: pub_month = pub_date['Month']
            except: pass
            try: pub_day = pub_date['Day']
            except: pass
            try: pages = article['Pagination'].get('StartPage', 'No pages')
            except: pass
            try: citation = f'{author_str}. {title}. {journal_title}. {pub_year}. {pub_month}. {pub_day};{journal_volume}({journal_issue}):{pages}. PMID: {pmid}.'
            except: citation = f'{title}. {journal_title}. PMID: {pmid}.'
            citation_arr.append(citation)
            abstract_arr.append(abstract_text)
            print(citation)
            print(abstract_text)
            print()
            print()
            print()
            create_folder(f'{query_folderpath}/yesterday/{pmid}')
            with open(f'{query_folderpath}/yesterday/{pmid}/abstract.txt', 'w') as f: f.write(abstract_text)
            row = [
                pmid,
                pub_year,
                pub_month,
                pub_day,
                journal_title,
                journal_volume,
                journal_issue,
                pages,
                title,
                author_str,
            ]
            time.sleep(random.randint(1, 3))
            util.csv_add_rows(f'{query_folderpath}/yesterday.csv', [row])
    else:
        print('no article found')
    print(f'number of abstracts: {len(abstract_arr)}')
    print(f'number of citations: {len(citation_arr)}')

def scrape_pubmed(query, year=None):
    query_slug = query.strip().lower().replace(' ', '-')
    pmid_list = get_ids(query, year, retmax=9999)
    print(pmid_list)
    citation_arr = []
    abstract_arr = []
    if pmid_list:
        create_folder(query_folderpath)
        for pmid_i, pmid in enumerate(pmid_list):
            print(f'{pmid_i}/{len(pmid_list)}')
            if os.path.exists(f'{query_folderpath}/master.csv'):
                pmid_csv = util.csv_get_rows_by_entity(f'{query_folderpath}/master.csv', pmid, col_num=0)
            else:
                row = [
                    'pmid',
                    'pub_year',
                    'pub_month',
                    'pub_day',
                    'journal_title',
                    'journal_volume',
                    'journal_issue',
                    'pages',
                    'title',
                    'author_st',
                ]
                util.csv_add_rows(f'{query_folderpath}/master.csv', [row])
                pmid_csv = []
            if pmid_csv != []: 
                continue
            abstract_text = ''
            title = ''
            authors = ''
            author_str = ''
            journal_title = ''
            journal_volume = ''
            journal_issue = ''
            pub_date = ''
            pub_year = ''
            pub_month = ''
            pub_day = ''
            pages = ''
            try: details = fetch_details(pmid)
            except: continue
            try: abstract_text = details['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            except: continue
            article = details['PubmedArticle'][0]['MedlineCitation']['Article']
            journal = article['Journal']
            pubmed_data = details['PubmedArticle'][0]['PubmedData']
            title = article.get('ArticleTitle', 'No title available')
            try: authors = article['AuthorList']
            except: pass
            try: author_str = ', '.join([f"{a['LastName']} {a['ForeName'][0]}" for a in authors])
            except: pass
            try: journal_title = journal.get('Title', 'No journal title available')
            except: pass
            try: journal_volume = journal['JournalIssue'].get('Volume', 'No volume')
            except: pass
            try: journal_issue = journal['JournalIssue'].get('Issue', 'No issue')
            except: pass
            try: pub_date = article.get('ArticleDate', [{'Year': 'No year', 'Month': 'No month', 'Day': 'No day'}])[0]
            except: pass
            try: pub_year = pub_date['Year']
            except: pass
            try: pub_month = pub_date['Month']
            except: pass
            try: pub_day = pub_date['Day']
            except: pass
            try: pages = article['Pagination'].get('StartPage', 'No pages')
            except: pass
            try: citation = f'{author_str}. {title}. {journal_title}. {pub_year}. {pub_month}. {pub_day};{journal_volume}({journal_issue}):{pages}. PMID: {pmid}.'
            except: citation = f'{title}. {journal_title}. PMID: {pmid}.'
            citation_arr.append(citation)
            abstract_arr.append(abstract_text)
            print(citation)
            print(abstract_text)
            print()
            print()
            print()
            create_folder(f'{query_folderpath}/{year}/{pmid}')
            with open(f'{query_folderpath}/{year}/{pmid}/abstract.txt', 'w') as f: f.write(abstract_text)
            row = [
                pmid,
                pub_year,
                pub_month,
                pub_day,
                journal_title,
                journal_volume,
                journal_issue,
                pages,
                title,
                author_str,
            ]
            util.csv_add_rows(f'{query_folderpath}/master.csv', [row])
            time.sleep(random.randint(1, 3))
    else:
        print('no article found')
    print(f'number of abstracts: {len(abstract_arr)}')
    print(f'number of citations: {len(citation_arr)}')

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

def gen_summaries_yesterday():
    rows = util.csv_get_rows(f'{query_folderpath}/yesterday.csv')
    for row in rows[1:]:
        pmid = row[0]
        pub_year = row[1]
        pub_month = row[2]
        pub_day = row[3]
        journal_title = row[4]
        journal_volume = row[5]
        journal_issue = row[6]
        pages = row[7]
        title = row[8]
        author_str = row[9]
        with open(f'{query_folderpath}/yesterday/{pmid}/abstract.txt') as f: abstract_text = f.read()
        summary_text = ''
        summary_text += f'ID: {pub_year}/{pub_month}/{pmid}\n'
        summary_text += f'DATE: {pub_year}/{pub_month}/{pub_day}\n'
        summary_text += f'JOURNAL: {journal_title}.{journal_volume}.{journal_issue}\n\n'
        summary_text += f'TITLE: {title}\n\n'
        prompt = f'''
            Summarize of the ABSTRACT below in 1 paragraph, in a way that is easy to understand.
            Focus the summary on ozone.
            Write very short and simple sentences.
            Don't write titles and subtitles.
            Don't reply with a list.
            Write the summary in Italian.
            ## ABSTRACT
            {abstract_text}
        '''
        reply = util_ai.gen_reply(prompt, model)
        summary_text += f'SUMMARY: {reply}\n\n'
        with open(f'{query_folderpath}/yesterday/{pmid}/summary.txt', 'w') as f: f.write(summary_text)
        print(summary_text)

def gen_folder_with_all_articles(foldername):
    for filename in os.listdir(query_folderpath):
        folderpath_in = f'{query_folderpath}/{filename}' 
        if not filename[0].isdigit(): continue
        if os.path.isdir(folderpath_in):
            studies_foldernames = os.listdir(f'{folderpath_in}')
            for study_foldername in studies_foldernames:
                abstract_filepath_in = f'{folderpath_in}/{study_foldername}/abstract.txt'
                abstract_folderpath_in = f'{folderpath_in}/{study_foldername}'
                print(abstract_filepath_in)
                abstract_folderpath_out = f'{query_folderpath}/{foldername}/{study_foldername}'
                abstract_filepath_out = f'{abstract_folderpath_out}/abstract.txt'
                create_folder(abstract_folderpath_out)
                with open(abstract_filepath_in) as f: content = f.read()
                with open(abstract_filepath_out, 'w') as f: f.write(content)

def get_latest_published_studies(num=5):
    master_rows = util.csv_get_rows(f'{query_folderpath}/master.csv')
    studies_rows = []
    for i in range(999):
        if len(studies_rows) >= 5: break
        date_curr = datetime.datetime.now() - datetime.timedelta(i)
        today_year = date_curr.year
        today_month = date_curr.month
        today_day = date_curr.day
        today_year = str(today_year)
        if today_month < 10: today_month = f'0{today_month}' 
        else: today_month = str(today_month)
        if today_day < 10: today_day = f'0{today_day}' 
        else: today_day = str(today_day)
        for row in master_rows:
            row_year = row[1]
            row_month = row[2]
            row_day = row[3]
            if row_year == '2024':
                print(today_year, today_month, today_day)
                print(row_year, row_month, row_day)
                print()
            if today_year == row_year and today_month == row_month and today_day == row_day: 
                studies_rows.append(row)
    for row in studies_rows:
        pmid = row[0] 
        folderpath_in = f'{query_folderpath}/news/todo/{pmid}'
        filepath_in = f'{query_folderpath}/news/todo/{pmid}/abstract.txt'
        folderpath_out = f'{query_folderpath}/news/latest/{pmid}'
        filepath_out = f'{query_folderpath}/news/latest/{pmid}/abstract.txt'
        create_folder(folderpath_out)
        with open(filepath_in) as f: content = f.read()
        with open(filepath_out, 'w') as f: f.write(content)
        print(row)

def move_to_bad_abstracts(pmid):
    latest_dir = f'{query_folderpath}/news/latest/{pmid}'
    source_dir = f'{query_folderpath}/news/todo/{pmid}'
    target_dir = f'{query_folderpath}/news/bad-abstract/{pmid}'
    try: shutil.move(source_dir, target_dir)
    except: print(f'[INFO] can\'t move {pmid} folder: folder not found')
    try: shutil.rmtree(latest_dir)
    except: print(f'[INFO] can\'t delete {pmid} folder: folder not found')

def summarize_article(pmid):
    source_filepath = f'{query_folderpath}/master.csv'
    rows = util.csv_get_rows_by_entity(source_filepath, str(pmid), col_num=0)
    if rows == []: return
    row = rows[0]
    with open(f'{query_folderpath}/news/latest/{pmid}/abstract.txt') as f: abstract_text = f.read()
    print(row)
    print(abstract_text)
    prompt = f'''
        Summarize in 1 paragraph the following ABSTRACT using the GUIDELINES below.
        ## GUIDELINES
        Focus the summary on ozone.
        Write the summary in an easy and simple way to understand.
        Write the summary in Italian.
        ## ABSTRACT
        {abstract_text}
    '''
    reply = util_ai.gen_reply(prompt)


def write_article(pmid):
    source_filepath = f'{query_folderpath}/master.csv'
    rows = util.csv_get_rows_by_entity(source_filepath, str(pmid), col_num=0)
    if rows == []: return
    row = rows[0]
    with open(f'{query_folderpath}/news/latest/{pmid}/abstract.txt') as f: abstract_text = f.read()
    prompt = f'''
        Write a detailed article about the ABSTRACT below and using the following GUIDELINES.
        ## GUIDELINES
        The article is 5 paragraphs.
        The paragraphs are very detailed and long.
        Write each paragraph in about 60-80 words.
        Title the paragraphs as following: Paragrafo 1, Paragrafo 2, Paragrafo 3, Paragrafo 4, Paragrafo 5.
        In paragraph 1, write the introduction.
        In paragraph 2, write the methods.
        In paragraph 3, write the results.
        In paragraph 4, write the discussions.
        In paragraph 5, write the conclusions.
        Don't repeat yourself. Don't repeat the same info more than once.
        Reply in Italian.
        ## ABSTRACT
        {abstract_text}
    '''
    reply = util_ai.gen_reply(prompt)
for year in years:
    scrape_pubmed_jsons(query, year)
    # time.sleep(600)
quit()
gen_folder_with_all_articles('all')
for year in years:
    scrape_pubmed(query, year)
    # time.sleep(600)
write_article(39088742)
# summarize_article(39088742)
move_to_bad_abstracts(38863229)
get_latest_published_studies(num=5)
scrape_pubmed_yesterday(query, yesterday=yesterday)
gen_summaries_yesterday()

    
