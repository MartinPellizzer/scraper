import textwrap

import numpy as np
from sentence_transformers import SentenceTransformer
from wikipediaapi import Wikipedia

import util_ai

wiki = Wikipedia('RAGBot/0.0', 'en')
doc = wiki.page('Hayao_Miyazaki').text
paragraphs = doc.split('\n\n')

for i, p in enumerate(paragraphs):
    wrapped_text = textwrap.fill(p, width=100)

    print('----------------------------------')
    print(wrapped_text)
    print('----------------------------------')

model = SentenceTransformer('Alibaba-NLP/gte-base-en-v1.5', trust_remote_code=True)
docs_embed = model.encode(paragraphs, normalize_embeddings=True)
print(docs_embed.shape)
print(docs_embed[0])

query = 'what was studio ghibli\'s first film?'
query_embed = model.encode(query, normalize_embeddings=True)
print(query_embed.shape)

similarities = np.dot(docs_embed, query_embed.T)
print(similarities.shape)
print(similarities)

top_3_idx = np.argsort(similarities, axis=0)[-3:][::-1].tolist()
print(top_3_idx)

most_similar_documents = [paragraphs[idx] for idx in top_3_idx]

context = ''
for i, p in enumerate(most_similar_documents):
    wrapped_text = textwrap.fill(p, width=100)

    print('----------------------------------')
    print(wrapped_text)
    print('----------------------------------')
    context += wrapped_text + '\n\n'

prompt = f'''
use the following CONTEXT to answer th QUESTION at the end.
if you don't know the answer, just say that you don't know, don't try to make up an answer.
CONTEXT: {context}
QUESTION: {query}
'''

reply = util_ai.gen_reply(prompt)

