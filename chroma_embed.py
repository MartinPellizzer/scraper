import os

import chromadb
from chromadb.utils import embedding_functions

from oliark import json_read
import util_ai

model = '/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf'

device = 'cpu'
device = 'cuda'

sentence_transformer_ef = embedding_functions.SentenceTransformerEmbeddingFunction(model_name='all-mpnet-base-v2', device=device)

chroma_client = chromadb.PersistentClient(path='ozonedb')
collection = chroma_client.get_or_create_collection(name='ozone', embedding_function=sentence_transformer_ef)

def embed():
    documents_folderpath = f'/home/ubuntu/vault/ozonogroup/studies/pubmed/ozone/json'
    documents_filenames = os.listdir(documents_folderpath)
    for i, document_filename in enumerate(documents_filenames):
        print(f'{i}/{len(documents_filenames)}')
        document_filepath = f'{documents_folderpath}/{document_filename}'
        data = json_read(document_filepath)
        try: abstract_text = data['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
        except: continue
        abstract_text = ' '.join(abstract_text).replace('  ', ' ')
        pmid = document_filename.split('.')[0]
        metadata = {'pmid': pmid}

        collection.add(
            documents=[abstract_text],
            metadatas=[metadata],
            ids=[pmid],
        )


query = input('ask >> ')

reply = util_ai.gen_reply('test', model)

results = collection.query(query_texts=[query], n_results=5)
context = results['documents'][0]
print(results)

prompt = f'''
    Answer the followin QUESTION using the following CONTEXT:
    ## QUESTION
    {query}
    ## CONTEXT
    {context}
'''
print(prompt)
reply = util_ai.gen_reply(prompt, model)

