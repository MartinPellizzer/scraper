import os
import numpy as np

import chromadb
from chromadb.utils import embedding_functions
from langchain.text_splitter import RecursiveCharacterTextSplitter, SentenceTransformersTokenTextSplitter
from sentence_transformers import CrossEncoder

import util
import util_ai


vault = '/home/ubuntu/vault'
model = '/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf'

proj = 'ozonogroup'
collection_name = 'ozone'
pmid_foldername = 'all'

pubmed_folderpath = f'{vault}/{proj}/studies/pubmed/{collection_name}/{pmid_foldername}'
master_folderpath = f'{vault}/{proj}/studies/pubmed/{collection_name}'

sentence_transformer_ef = embedding_functions.SentenceTransformerEmbeddingFunction(model_name='all-mpnet-base-v2')

chroma_client = chromadb.PersistentClient(path='database')
collection = chroma_client.get_or_create_collection(name=collection_name, embedding_function=sentence_transformer_ef)

def add_embeddings():
    studies_folderpath = pubmed_folderpath
    studies_foldernames = os.listdir(studies_folderpath)
    for i, study_foldername in enumerate(studies_foldernames):
        if not os.path.isdir(f'{studies_folderpath}/{study_foldername}'): continue
        print(f'{i}/{len(studies_foldernames)} - {study_foldername}')
        study_abstract_filepath = f'{studies_folderpath}/{study_foldername}/abstract.txt'
        with open(study_abstract_filepath) as f: abstract_text = f.read()
        rows = util.csv_get_rows_by_entity(f'{master_folderpath}/master.csv', study_foldername, col_num=0)
        if rows == []: continue
        row = rows[0]
        pmid = row[0]
        year = row[1]
        month = row[2]
        day = row[3]
        journal_title = row[4]
        journal_volume = row[5]
        journal_issue = row[6]
        title = row[7]
    
        collection.add(
            documents=[abstract_text],
            metadatas=[
                {'pmid': pmid,
                'year': year,
                'month': month,
                'day': day,
                'journal_title': journal_title,
                'journal_volume': journal_volume,
                'journal_issue': journal_issue,
                'title': title}
            ],
            ids=[study_foldername]
        )

# add_embeddings()
# quit()

def gen_aux_questions(query):
    # [augment] generate auxiliary questions
    prompt = f'''
        You are a helpful scientific paper writer. 
        Suggest up to 5 additional related questions to the provided QUESTION, to help find the information needed for the provided question.
        Suggest only short questions without compound sentences.
        Suggest a variety of questions that cover different aspects of the topic.
        Make sure they are completed questions, and that they are related to the original question.
        Output one question per line. Do not number the questions.
        QUESTION: {query}
    '''.strip()
    reply = util_ai.gen_reply(prompt, model)
    related_questions = [related_question for related_question in reply.split('\n') if related_question.strip() != ''] 
    return related_questions


def gen_ex_answers(questions):
    # [augment] generate example answers
    queries = []
    for question in questions:
        prompt = f'''
        You are a helpful scientific paper writer. Provide an example answer to the following QUESTION, that may be found in a scientific paper abstract. Reply only with the abstract. Don't include the title.
        QUESTION: {question}
        '''.strip()
        reply = util_ai.gen_reply(prompt, model)
        queries.append(f'{query} {question} {reply}')
    for question in questions:
        print('###############################################')
        print(question)
        print('###############################################')
    return queries

def retrieve_docs(queries):
    # [get] get docs from vectordb
    n_results = 5
    n_queries = len(queries)
    results = collection.query(query_texts=queries, n_results=n_results, include=['documents', 'metadatas'])
    print(results['ids'])
    unique_results = []
    for i in range(n_queries):
        for j in range(n_results):
            dict_tmp = {}
            dict_tmp['metadata'] = results['metadatas'][i][j]
            dict_tmp['document'] = results['documents'][i][j]
            found = False
            for unique_result in unique_results:
                if dict_tmp['metadata']['pmid'] == unique_result['metadata']['pmid']:
                    found = True
                    break
            if not found:
                unique_results.append(dict_tmp)
    for document in unique_results:
        print(document)
        print()
    print(len(unique_results))
    return unique_results

def order_retrieved_docs(unique_results):
    # [order] order retrieved documents by most relevant
    cross_encoder = CrossEncoder('cross-encoder/ms-marco-MiniLM-L-6-v2')
    unique_documents = [result['document'] for result in unique_results]
    for document in unique_documents:
        print(document)
        print()
    pairs = [[query, doc] for doc in unique_documents]
    scores = cross_encoder.predict(pairs)
    for score in scores:
        print(score)
    ordered_results = []
    for o in np.argsort(scores)[::-1]:
        print(o)
        ordered_results.append(int(o))
    print(ordered_results)
    return ordered_results

question = ''
while True:
    query = input('ask grim >>> ')
    if query != '': question = query
    if question == '': continue
    related_questions = gen_aux_questions(question)
    queries = gen_ex_answers(related_questions)
    unique_results = retrieve_docs(queries)
    ordered_results = order_retrieved_docs(unique_results)

    print('##############################################################')
    print('##############################################################')
    print('##############################################################')
    best_results = []
    for i in range(5):
        best_results.append(unique_results[ordered_results[i]])
        print(unique_results[ordered_results[i]])
        print()
    print('##############################################################')
    print('##############################################################')
    print('##############################################################')
    
    context = ''
    for result in best_results:
        document = result['document']
        context += f'{document}\n\n'
    prompt = f'''
        Answer the followin QUESTION using the following CONTEXT:
        ## QUESTION
        {question}
        ## CONTEXT
        {context}
    '''
    print(prompt)
    reply = util_ai.gen_reply(prompt, model)
    if 'not enough context' not in reply.strip().lower():
        done = True
        
    print('##############################################################')
    print(reply)
    print('##############################################################')

