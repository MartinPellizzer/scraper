import os 
import shutil
import datetime
import requests
import base64

from PIL import Image, ImageFont, ImageDraw

import util
import util_ai

model = '/Meta-Llama-3.1-8B-Instruct-Q4_K_M.gguf'

proj = 'ozonogroup'
query = f'ozone'.strip().lower()
query_slug = query.replace(' ', '-')

vault = '/home/ubuntu/vault'
pubmed_folderpath = f'{vault}/{proj}/studies/pubmed'
query_folderpath = f'{pubmed_folderpath}/{query_slug}'
master_filepath = f'{query_folderpath}/master.csv'
news_folderpath = f'{vault}/{proj}/news'

def get_latest_published_study_json():
    for i in range(99999):
        date_curr = datetime.datetime.now() - datetime.timedelta(i)
        today_year = date_curr.year
        today_month = date_curr.month
        today_day = date_curr.day
        today_year = str(today_year)
        if today_month < 10: today_month = f'0{today_month}' 
        else: today_month = str(today_month)
        if today_day < 10: today_day = f'0{today_day}' 
        else: today_day = str(today_day)
        jsons_filenames_todo = [f for f in os.listdir(f'{query_folderpath}/json')]
        jsons_filenames_done = [f for f in os.listdir(f'{news_folderpath}/done')]
        jsons_filenames_rejected = [f for f in os.listdir(f'{news_folderpath}/rejected')]
        for json_filename_todo in jsons_filenames_todo:
            if json_filename_todo in jsons_filenames_done: continue
            if json_filename_todo in jsons_filenames_rejected: continue
            json_filepath = f'{query_folderpath}/json/{json_filename_todo}'
            data = util.json_read(json_filepath)
            try: citation = data['PubmedArticle'][0]['MedlineCitation']
            except: continue
            try: pmid = citation['PMID']
            except: continue
            try: article = data['PubmedArticle'][0]['MedlineCitation']['Article']
            except: continue
            try: date = article['ArticleDate'][0]
            except: continue
            try: year = date['Year']
            except: continue
            try: month = date['Month']
            except: continue
            try: day = date['Day']
            except: continue
            try: title = article['ArticleTitle']
            except: continue
            try: abstract_text = article['Abstract']['AbstractText']
            except: continue
            try: journal_title = article['Journal']['Title']
            except: journal_title = ''
            try: journal_volume = article['Journal']['JournalIssue']['Volume']
            except: journal_volume
            try: journal_issue = article['Journal']['JournalIssue']['Issue']
            except: journal_issue
            abstract_text = ' '.join(abstract_text)
            if today_year == year and today_month == month and today_day == day: 
                if False:
                    print(f'[PMID] {pmid}')
                    print(f'[DATE] {year}/{month}/{day}')
                    print(f'[JOURNAL] {journal_title}.{journal_volume}.{journal_issue}')
                    print(f'[TITLE] {title}')
                    print(f'[ABSTRACT] {abstract_text}')
                return {
                    'pmid': pmid,
                    'year': year,
                    'month': month,
                    'day': day,
                    'journal_title': journal_title,
                    'journal_volume': journal_volume,
                    'journal_issue': journal_issue,
                    'study_title': title,
                    'abstract_text': abstract_text,
                }

def write_article_json(data):
    pmid = data['pmid']
    year = data['year']
    month = data['month']
    day = data['day']
    journal_title = data['journal_title']
    journal_volume = data['journal_volume']
    journal_issue = data['journal_issue']
    study_title = data['study_title']
    abstract_text = data['abstract_text']
    prompt = f'''test'''
    reply = util_ai.gen_reply(prompt)
    print()
    print(pmid)
    print(f'{year}/{month}/{day}')
    print()
    print(study_title)
    print()
    print(abstract_text)
    print()
    prompt = f'''
        Translate in Italian the following ABSTRACT.
        ## ABSTRACT
        {study_title}
        {abstract_text}
    '''
    reply = util_ai.gen_reply(prompt)
    for i in range(3):
        if i >= 3: return
        is_abstract_good = input('\ngood abstract? [y/n] >>> ')
        if is_abstract_good == 'y': break
        if is_abstract_good == 'n': 
            json_filepath_source = f'{pubmed_folderpath}/{query_slug}/json/{pmid}.json'
            json_filepath_target = f'{news_folderpath}/rejected/{pmid}.json'
            shutil.copy2(json_filepath_source, json_filepath_target)
            return
    
    json_filepath = f'{news_folderpath}/curr/{pmid}.json'
    util.create_folder_for_filepath(json_filepath)
    util.json_generate_if_not_exists(json_filepath)
    data = util.json_read(json_filepath)
    data['id'] = pmid
    data['year'] = year
    data['month'] = month
    data['day'] = day
    util.json_write(json_filepath, data)

    key = 'title'
    if key not in data: data[key] = ''
    if data[key] == '':
        run = True
        while run:
            prompt = f'''
                Write a numbered list of 10 titles for an article about the ABSTRACT below.
                Use the following GUIDELINES to write the titles.
                ## GUIDELINES
                Each headline must be unique, novel, and different from the others.
                Each headline must be 3 to 16 words long.
                Each headline must focus on ozone.
                Reply in Italian.
                ## ABSTRACT
                {study_title}
                {abstract_text}
            '''
            reply = util_ai.gen_reply(prompt)
            lines = []
            for line in reply.split('\n'):
                line = line.strip()
                if line == '': continue
                if not line[0].isdigit(): continue
                if '. ' not in line: continue
                line = '. '.join(line.split('. ')[1:])
                line = line.strip()
                if line == '': continue
                line = line.replace('.', '') 
                line = line.replace('"', '') 
                line = line.replace('**', '') 
                lines.append(line)
            if len(lines) == 10:
                print('*******************************************')
                print('*******************************************')
                for i in range(len(lines)):
                    print(f'{i}: {lines[i]}')
                selected = input('select title num >>> ')
                if selected.strip() == '': continue
                selected = int(selected)
                data['title'] = lines[selected]
                util.json_write(json_filepath, data)
                run = False
    key = 'slug'
    if key not in data: data[key] = ''
    if data[key] == '':
        run = True
        while run:
            title = data['title']
            prompt = f'''
                Rewrite the following title by removing all stop words, determinative articles, and other unnecessary words.
                Keep only the most important words.
                Keep only the keywords.
                Keep the spaces.
                Write as few words as possible.
                The end results should look like a "slug" you find in a url.
                Reply only with the slug.
                ## TITLE
                {title}
            '''
            reply = util_ai.gen_reply(prompt)
            slug = reply.replace(' ', '-').lower().strip()
            print('************************************************')
            print(slug)
            print('************************************************')
            slug_ok = input('press y to confirm slug >>> ')
            if slug_ok == 'y':
                data['slug'] = slug
                util.json_write(json_filepath, data)
                run = False
    key = 'category'
    if key not in data: data[key] = ''
    if data[key] == '':
        categories = ['sanificazione', 'medicina', 'ambiente', 'lavorazione', 'tecnologia', 'chimica']
        for i, category in enumerate(categories):
            print(f'{i}. {category}')
        selected = int(input('select title num >>> '))
        data[key] = categories[selected]
        util.json_write(json_filepath, data)
    key = 'body'
    if key not in data: data[key] = []
    if data[key] == []:
        run = True
        while run:
            prompt = f'''
                Write a detailed article discussing the scientific study in the ABSTRACT below and using the following GUIDELINES.
                ## GUIDELINES
                The article is 5 paragraphs.
                The paragraphs are very detailed and long.
                Write each paragraph in about 60-80 words.
                Title the paragraphs as following: Paragrafo 1, Paragrafo 2, Paragrafo 3, Paragrafo 4, Paragrafo 5.
                In paragraph 1, write the introduction and iclude the fact that the scientific study discussed was publisced by the following publication: "{journal_title}".
                In paragraph 2, write the methods.
                In paragraph 3, write the results.
                In paragraph 4, write the discussions.
                In paragraph 5, write the conclusions.
                Write the paragraphs in an easy and simple to understand way.
                Clarify and expand on concepts that are not easily understandable by most people.
                Include the name of the publication in the article: {journal_title}. 
                Include info about ozone from the abstract.
                Don't repeat yourself.
                Don't include the title of the study in the abstract.
                Don't include refereces to links, e.g. (1), etc...
                Reply in Italian.
                ## ABSTRACT
                {data['title']}
                {abstract_text}
            '''
            reply = util_ai.gen_reply(prompt)
            paragraphs = []
            for line in reply.split('\n'):
                line = line.strip()
                if line == '': continue
                if 'paragrafo' in line.lower(): continue
                if '**' in line.lower(): continue
                paragraphs.append(line)
            if len(paragraphs) == 5:
                print('*******************************************')
                for paragraph in paragraphs:
                    print(paragraph)
                    print()
                print('*******************************************')
                confirmed = input('press y to confirm body >>> ')
                if confirmed == 'y':
                    data[key] = paragraphs
                    util.json_write(json_filepath, data)
                    run = False
    # img
    run = True
    if not os.path.exists(f'{news_folderpath}/curr/{pmid}-raw.jpg'):
        while run:
            print(data['body'])
            prompt = input('insert image promt >>> ')
            if prompt.strip() == '': continue
            payload = {
                "prompt": prompt,
                "width": 1024,
                "height": 1024,
                "steps": 25,
                "cfg_scale": 6,
                "denoising_strength": 0.7,
                "sampler_name": "DPM++ 2M",
                "scheduler": "Karras",
                "seed": -1,
                "batch_size": 1
            }
            response = requests.post(url='http://127.0.0.1:7860/sdapi/v1/txt2img', json=payload)
            r = response.json()
            print(r)
            export_filepath = f'{news_folderpath}/tmp/{pmid}.jpg'
            with open(export_filepath, 'wb') as f:
                f.write(base64.b64decode(r['images'][0]))
            img = Image.open(export_filepath)
            img.show()
            confirmed = input('confirm image? [y] >>> ')
            if confirmed == 'y':
                export_filepath = f'{news_folderpath}/curr/{pmid}-raw.jpg'
                with open(export_filepath, 'wb') as f:
                    f.write(base64.b64decode(r['images'][0]))
                img = Image.open(f'{news_folderpath}/curr/{pmid}-raw.jpg')
                img.thumbnail((768, 768), Image.Resampling.LANCZOS)
                draw = ImageDraw.Draw(img)
                font = ImageFont.truetype('fonts/Lato/Lato-Regular.ttf', 16)
                text = 'ozonogroup'
                _, _, text_w, text_h = font.getbbox(text)
                padding = int(768 * 0.05)
                draw.text((768 - text_w - padding, 768 - text_h - padding), text, '#ffffff', font=font)
                img.save(f'{news_folderpath}/curr/{pmid}-final.jpg', optimize=True, quality=70)
                run = False
            if confirmed == 's':
                run = False

    if False:
        if os.path.exists(f'{news_folderpath}/curr/{pmid}-raw.jpg'):
            if not os.path.exists(f'{news_folderpath}/curr/{pmid}-final.jpg'):
                confirm = input('process img? [y] >>> ')
                if confirm == 'y':
                    img = Image.open(f'{news_folderpath}/curr/{pmid}-raw.jpg')
                    img.thumbnail((768, 768), Image.Resampling.LANCZOS)
                    draw = ImageDraw.Draw(img)
                    font = ImageFont.truetype('fonts/Lato/Lato-Regular.ttf', 16)
                    text = 'ozonogroup'
                    _, _, text_w, text_h = font.getbbox(text)
                    padding = int(768 * 0.05)
                    draw.text((768 - text_w - padding, 768 - text_h - padding), text, '#ffffff', font=font)
                    img.save(f'{news_folderpath}/curr/{pmid}-test.jpg', optimize=True, quality=70)

    confirmed = input('confirm article? [y] >>> ')
    if confirmed == 'y':
        json_filepath_done = f'{news_folderpath}/done/{pmid}.json'
        image_filepath_source = f'{news_folderpath}/curr/{pmid}-final.jpg'
        image_filepath_target = f'{news_folderpath}/images/{pmid}.jpg'
        shutil.move(json_filepath, json_filepath_done)
        shutil.move(image_filepath_source, image_filepath_target)

data = get_latest_published_study_json()
write_article_json(data)
