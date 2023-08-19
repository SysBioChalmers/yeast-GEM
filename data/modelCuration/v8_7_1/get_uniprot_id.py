import requests
from bs4 import BeautifulSoup
import time
from tqdm import tqdm, trange
import re
import pandas as pd
import time

def crawler_uniprot(gene):
    headers = {'user-agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36"}
    url = r'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&query=%28{}%29'.format(gene)
    response = requests.request("GET", url, headers=headers)
    rawdata = BeautifulSoup(response.text, 'html.parser')
    try:
        recom0 = re.compile(r'"primaryAccession".*S288c')
        id0 = recom0.findall(str(rawdata))
        recom = re.compile(r'"primaryAccession":"[0-9a-zA-Z]*"')
        id = recom.findall(id0[0])
        uniprotid = id[0].split('"')[3]
    except:
        uniprotid = ""
    return uniprotid

sgd = pd.read_table('../../databases/SGDgeneNames.tsv')
sgd.index = sgd.loc[:, 'Systematic_name']
genename = list(sgd.loc[:, 'Systematic_name'])
c = 0
for g in genename:
    try:
        time.sleep(0.5)
        print(c/len(genename))
        print(g)
        sgd.loc[g, 'Uniprot_id'] = crawler_uniprot(g)
        print(sgd.loc[g, 'Uniprot_id'])
        c += 1
    except requests.exceptions.ConnectionError:
        time.sleep(5)
        print(c/len(genename))
        print(g)
        sgd.loc[g, 'Uniprot_id'] = crawler_uniprot(g)
        print(sgd.loc[g, 'Uniprot_id'])
        c += 1
sgd.to_csv('./SGD_with_Uniprot.csv', index=False)


