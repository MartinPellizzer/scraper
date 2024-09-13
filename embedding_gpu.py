import os
import time
import json

from tqdm.auto import tqdm

from sentence_transformers import SentenceTransformer

proj = 'ozonogroup'
vault = '/home/ubuntu/vault'

proj_folderpath = f'{vault}/{proj}'
studies_folderpath = f'{proj_folderpath}/studies'
pubmed_folderpath = f'{studies_folderpath}/pubmed'
ozone_folderpath = f'{pubmed_folderpath}/ozone'
json_folderpath = f'{ozone_folderpath}/json'

jsons_filepaths = [f'{json_folderpath}/{filename}' for filename in os.listdir(json_folderpath)]

abstracts = []
for json_filepath in jsons_filepaths:
    print(json_filepath)
    with open(json_filepath) as f: data = json.load(f)
    try: article = data['PubmedArticle'][0]['MedlineCitation']['Article']
    except: continue
    try: abstract_text = article['Abstract']['AbstractText']
    except: continue
    abstract_text = ' '.join(abstract_text)
    abstracts.append(abstract_text)

embedding_model = SentenceTransformer(model_name_or_path='all-mpnet-base-v2', device='cpu')

sentences = abstracts[:3]

if False:
    start_time = time.time()

    embedding_model.to('cuda')
    for item in tqdm(abstracts):
        embedding = embedding_model.encode(item)

    end_time = time.time()
    delta_time = end_time - start_time

    print(f'--- {delta_time} ---')

if True:
    start_time = time.time()

    embedding_model.to('cuda')
    chunk_embeddings = embedding_model.encode(abstracts[:1000], batch_size=4, convert_to_tensor=True)

    end_time = time.time()
    delta_time = end_time - start_time

    print(f'--- {delta_time} ---')

'''
embeddings_dict = dict(zip(sentences, embeddings))
for sentence, embedding in embeddings_dict.items():
    print(sentence)
    print(embedding)
    print()
'''
