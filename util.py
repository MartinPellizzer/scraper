import os
import csv
import json

def create_folder_for_filepath(filepath):
    chunks = filepath.split('/')
    chunk_curr = ''
    for chunk in chunks[:-1]:
        chunk_curr += chunk + '/'
        try: os.makedirs(chunk_curr)
        except: pass


###################################
# CSV
###################################

def csv_get_rows(filepath, delimiter='\\'):
    rows = []
    with open(filepath, encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for i, line in enumerate(reader):
            rows.append(line)
    return rows


def csv_add_rows(filepath, rows, delimiter='\\'):
    with open(filepath, 'a', encoding='utf-8', errors='ignore', newline='') as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerows(rows)


def csv_set_rows(filepath, rows, delimiter='\\'):
    with open(filepath, 'w', encoding='utf-8', errors='ignore', newline='') as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerows(rows)
        

def csv_get_rows_by_entity(filepath, entity, delimiter='\\', col_num=0):
    rows = []
    with open(filepath, encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for i, line in enumerate(reader):
            if line == []: continue
            if line[col_num].lower().strip() == entity.lower().strip():
                rows.append(line)
    return rows


def folder_create(path):
    if not os.path.exists(path): os.makedirs(path)


def csv_get_header_dict(rows):
    cols = {}
    for i, val in enumerate(rows[0]):
        cols[val] = i
    return cols


def csv_get_cols(rows):
    cols = {}
    for i, val in enumerate(rows[0]):
        cols[val] = i
    return cols


def csv_get_rows_filtered(filepath, col_num, col_val, delimiter='\\'):
    rows = []
    with open(filepath, encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for i, line in enumerate(reader):
            if line == []: continue
            if line[col_num].lower().strip() == col_val.lower().strip():
                rows.append(line)
    return rows




###################################
# FILE
###################################

def file_read(filepath):
    with open(filepath, 'a', encoding='utf-8') as f: pass
    with open(filepath, 'r', encoding='utf-8') as f: 
        text = f.read()
    return text


def file_append(filepath, text):
    with open(filepath, 'a', encoding='utf-8') as f: 
        f.write(text)


def file_write(filepath, text):
    create_folder_for_filepath(filepath)
    with open(filepath, 'w', encoding='utf-8') as f: f.write(text)





###################################
# JSON
###################################

def json_generate_if_not_exists(filepath):
    if not os.path.exists(filepath):
        file_append(filepath, '')

    if file_read(filepath).strip() == '':
        file_append(filepath, '{}')


def json_append(filepath, data):
    with open(filepath, 'a', encoding='utf-8') as f:
        json.dump(data, f)


def json_read(filepath):
    with open(filepath, 'r', encoding='utf-8') as f: 
        return json.load(f)


def json_write(filepath, data):
    with open(filepath, 'w', encoding='utf-8') as f:
        json.dump(data, f)





