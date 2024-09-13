import os
import time
import textwrap

from Bio import Entrez
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sentence_transformers import SentenceTransformer

import util
import util_ai

sort_by = 'pub_date'
sort_by = 'relevance'

vault = '/home/ubuntu/vault'

model = 'Mistral-Nemo-Instruct-2407.Q4_0.gguf'
model = '/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf'

query = 'ozone'
year = '2023'

export_folder = f'{vault}/ozonogroup/studies/pubmed/{year}'
export_folder_ids = f'{vault}/ozonogroup/studies/pubmed'
try: os.makedirs(f'{vault}/ozonogroup/studies')
except: pass
try: os.makedirs(f'{vault}/ozonogroup/studies/pubmed')
except: pass
try: os.makedirs(f'{vault}/ozonogroup/studies/pubmed/{year}')
except: pass

Entrez.email = 'leenrandell@gmail.com'

def get_ids(query, year):
    handle = Entrez.esearch(db='pubmed', term=query, retmax=9999, sort=sort_by, datetype='edat', mindate=year, maxdate=year)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def fetch_details(pmid):
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    records = Entrez.read(handle)
    handle.close()
    return records

def scrape_pubmed():
    pmid_list = get_ids(query, year)
    citation_arr = []
    abstract_arr = []
    if pmid_list:
        ids = '\n'.join(pmid_list)
        with open(f'{export_folder_ids}/{year}-ids.txt', 'w') as f: f.write(ids)
        for pmid_i, pmid in enumerate(pmid_list):
            print(f'{pmid_i}/{len(pmid_list)}')
            pmid_csv = util.csv_get_rows_by_entity(f'{vault}/ozonogroup/studies/pubmed/master.csv', pmid, col_num=0)
            if pmid_csv != []: continue
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
            details = fetch_details(pmid)
            try:
                abstract_text = details['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            except:
                continue
            article = details['PubmedArticle'][0]['MedlineCitation']['Article']
            journal = article['Journal']
            pubmed_data = details['PubmedArticle'][0]['PubmedData']
            title = article.get('ArticleTitle', 'No title available')
            try:
                authors = article['AuthorList']
            except: pass
            try:
                author_str = ', '.join([f"{a['LastName']} {a['ForeName'][0]}" for a in authors])
            except: pass
            try:
                journal_title = journal.get('Title', 'No journal title available')
            except: pass
            try:
                journal_volume = journal['JournalIssue'].get('Volume', 'No volume')
            except: pass
            try:
                journal_issue = journal['JournalIssue'].get('Issue', 'No issue')
            except: pass
            try:
                pub_date = article.get('ArticleDate', [{'Year': 'No year', 'Month': 'No month', 'Day': 'No day'}])[0]
            except: pass
            try:
                pub_year = pub_date['Year']
            except: pass
            try:
                pub_month = pub_date['Month']
            except: pass
            try:
                pub_day = pub_date['Day']
            except: pass
            try:
                pages = article['Pagination'].get('StartPage', 'No pages')
            except: pass
            try:
                citation = f'{author_str}. {title}. {journal_title}. {pub_year}. {pub_month}. {pub_day};{journal_volume}({journal_issue}):{pages}. PMID: {pmid}.'
            except:
                citation = f'{title}. {journal_title}. PMID: {pmid}.'
            citation_arr.append(citation)
            abstract_arr.append(abstract_text)

            print(citation)
            print(abstract_text)
            print()
            print()
            print()
            try: os.makedirs(f'{vault}/ozonogroup/studies/pubmed/{year}/{pmid}')
            except: pass
            with open(f'{export_folder}/{pmid}/abstract.txt', 'w') as f: f.write(abstract_text)
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
                author_str
            ]
            util.csv_add_rows(f'{vault}/ozonogroup/studies/pubmed/master.csv', [row])
            time.sleep(1)
    else:
        print('no article found')

    print(f'number of abstracts: {len(abstract_arr)}')
    print(f'number of citations: {len(citation_arr)}')
    abstracts_text = '\n\n'.join(abstract_arr)
    
    return citation_arr, abstract_arr

scrape_pubmed()
quit()

def rag():
    # citation_arr, abstract_arr = scrape_pubmed()
    abstract_arr = []
    studies_folderpath = f'{vault}/studies/terrawhisper/pubmed'
    studies_filenames = os.listdir(studies_folderpath)
    for study_filename in studies_filenames:
        study_filepath = f'{studies_folderpath}/{study_filename}'
        with open(study_filepath) as f: abstract_text = f.read()
        abstract_arr.append(abstract_text)

    model = SentenceTransformer('Alibaba-NLP/gte-base-en-v1.5', trust_remote_code=True)
    docs_embed = model.encode(abstract_arr, normalize_embeddings=True)
    print(docs_embed.shape)
    print(docs_embed[0])

    query_embed = model.encode(query, normalize_embeddings=True)
    print(query_embed.shape)

    similarities = np.dot(docs_embed, query_embed.T)
    print(similarities.shape)
    print(similarities)

    top_3_idx = np.argsort(similarities, axis=0)[-3:][::-1].tolist()
    print(top_3_idx)

    most_similar_documents = [abstract_arr[idx] for idx in top_3_idx]
    abstracts_text = '\n\n'.join(most_similar_documents)

    context = ''
    for i, p in enumerate(most_similar_documents):
        wrapped_text = textwrap.fill(p, width=100)

        print('----------------------------------')
        print(wrapped_text)
        print('----------------------------------')
        context += wrapped_text + '\n\n'

    prompt = f'''
    write a summary using the data from the ABSTRACTS below.
    also, follow the GUIDELINES below.
    GUIDELINES:
    reply only with the paragraph
    ABSTRACTS: 
    {abstracts_text}
    '''

    reply = util_ai.gen_reply(prompt)

def cluster():
    abstract_arr = []
    studies_folderpath = export_folder
    studies_filenames = os.listdir(studies_folderpath)
    for study_filename in studies_filenames[:200]:
        study_filepath = f'{studies_folderpath}/{study_filename}'
        with open(study_filepath) as f: abstract_text = f.read()
        abstract_arr.append(abstract_text)

    model = SentenceTransformer('Alibaba-NLP/gte-base-en-v1.5', trust_remote_code=True)
    docs_embed = model.encode(abstract_arr, normalize_embeddings=True)
    print(docs_embed.shape)
    print(docs_embed[0])
    numeric_vectors = docs_embed

    wcss = []
    for i in range(1, len(numeric_vectors)):
        kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
        kmeans.fit(numeric_vectors)
        wcss.append(kmeans.inertia_)

    plt.figure(figsize=(10,5))
    plt.plot(range(1, len(numeric_vectors)), wcss)
    plt.title('elbow method')
    plt.xlabel('num of clusters')
    plt.ylabel('wcss')
    plt.show()
    
cluster()
quit()


rag()










quit()
with open('abstracts-test.txt') as f: abstracts_text = f.read()

'''
model = SentenceTransformer('Alibaba-NLP/gte-base-en-v1.5', trust_remote_code=True)
docs_embed = model.encode(abstract_arr, normalize_embeddings=True)
print(docs_embed.shape)
print(docs_embed[0])
'''
def gen_herbs():
    prompt = f'''
    give me a numbered list of medicinal herbs that helps with cough.
    follow the GUIDELINES below.
    get the herbs from the ABSTRACTS below.
    GUIDELINES:
    - reply with only the names of the herbs mentioned in the ABSTRACTS, don't make up names
    ABSTRACTS:
    {abstracts_text}
    '''
    reply = util_ai.gen_reply(prompt, model)
