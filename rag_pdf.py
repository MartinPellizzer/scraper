import random

import fitz
from tqdm.auto import tqdm
import pandas as pd
from spacy.lang.en import English

pdf_path = 'ozone-food-processing.pdf'

def pretty(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            print('\t' * (indent+1) + str(value))
            

def text_formatter(text):
    cleaned_text = text.replace('\n', ' ').strip()
    return cleaned_text

def open_and_read_pdf(pdf_path):
    doc = fitz.open(pdf_path)
    pages_and_texts = []
    for page_number, page in tqdm(enumerate(doc)):
        text = page.get_text()
        text = text_formatter(text=text)
        pages_and_texts.append({
            'page_number': page_number+1,
            'page_char_count': len(text),
            'page_word_count': len(text.split(' ')),
            'page_sentence_count_raw': len(text.split('. ')),
            'page_token_count': len(text) / 4,
            'text': text,
        })
    return pages_and_texts

pages_and_texts = open_and_read_pdf(pdf_path=pdf_path)

def preview_pages_and_texts(pages_and_texts):
    print(random.sample(pages_and_texts, k=3))
    df = pd.DataFrame(pages_and_texts)
    print(df.head())
    print(df.describe().round(2))

nlp = English()
nlp.add_pipe('sentencizer')
for item in tqdm(pages_and_texts):
    item['sentences'] = list(nlp(item['text']).sents)
    item['sentences'] = [str(sentence) for sentence in item['sentences']]
    item['page_sentence_count_spacy'] = len(item['sentences'])

pretty(random.sample(pages_and_texts, k=1)[0])

num_sentence_chunk_size = 10

def split_list(input_list, slice_size=num_sentence_chunk_size):
    return [input_list[i:i+slice_size] for i in range(0, len(input_list ), slice_size)]

test_list = list(range(25))
print(split_list(test_list))

for item in tqdm(pages_and_texts):
    item['sentence_chunks'] = split_list(input_list=item['sentences'], slice_size=num_sentence_chunk_size)
    item['num_chunks'] = len(item['sentence_chunks'])

# pretty(random.sample(pages_and_texts, k=1)[0])
pretty(pages_and_texts[21-1])
df = pd.DataFrame(pages_and_texts)
print(df.describe().round(2))


## split each chunk into its own item

import re
pages_and_chunks = []
for item in tqdm(pages_and_texts):
    for sentence_chunk in item['sentence_chunks']:
        chunk_dict = {}
        chunk_dict['page_number'] = item['page_number']
        joined_sentence_chunk = ''.join(sentence_chunk).replace('  ', ' ').strip()
        joined_sentence_chunk = re.sub(r'\.([A-Z])', r'. \1', joined_sentence_chunk) 
        chunk_dict['sentence_chunk'] = joined_sentence_chunk
        chunk_dict['chunk_char_count'] = len(joined_sentence_chunk)
        chunk_dict['chunk_word_count'] = len([word for word in joined_sentence_chunk.split(' ')])
        chunk_dict['chunk_token_count'] = len(joined_sentence_chunk) / 4
        pages_and_chunks.append(chunk_dict)

print(len(pages_and_chunks))
        
df = pd.DataFrame(pages_and_chunks)
print(df.describe().round(2))

## filter out chunks with too few token
min_token_length = 30
for row in df[df['chunk_token_count'] <= min_token_length].sample(5).iterrows():
    print(f'{row[1]["sentence_chunk"]}')

pages_and_chunks_over_min_token_len = df[df['chunk_token_count'] > min_token_length].to_dict(orient='records')
print(pages_and_chunks_over_min_token_len[:2])
