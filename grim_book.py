
import fitz
from tqdm.auto import tqdm
from langchain.text_splitter import RecursiveCharacterTextSplitter

pdf_path = 'ozone-food-processing.pdf'

book_text = ''
doc = fitz.open(pdf_path)
for page_number, page in tqdm(enumerate(doc)):
    text = page.get_text().strip()
    book_text += f'{text} '

'''
text_splitter = RecursiveCharacterTextSplitter(chunk_size=384*4, chunk_overlap=0)
chunks = text_splitter.create_documents([book_text])

for chunk in chunks:
    print(chunk)

print(type(chunks))
'''

'''
import os
from unstructured.partition.pdf import partition_pdf
from unstructured.staging.base import element_to_json

elements = partition_pdf(
    filename=pdf_path,
    strategy='hi_res',
    infer_table_structure=True,
    model_name='yolox',
)
'''

## semantic chunking
import re

single_sentences_list = re.split(f'(?<=[.?!])\s+', book_text)
print(len(single_sentences_list))

sentences = [{'sentence': x, 'index': i} for i, x in enumerate(single_sentences_list)]

def combine_sentences(sentences, buffer_size=1):
    for i in range(len(sentences)):
        combined_sentence = ''
        for j in range(i - buffer_size, i):
            if j >= 0:
                combined_sentence += sentences[j]['sentence'] + ' '

        combined_sentence += sentences[i]['sentence']
        for j in range(i + 1, i + 1 + buffer_size):
            if j < len(sentences):
                combined_sentence += ' ' + sentences[j]['sentence']

        sentences[i]['combined_sentence'] = combined_sentence

    return sentences

sentences = combine_sentences(sentences)
for sentence in sentences[:3]:
    print(sentence)
