import requests
from bs4 import BeautifulSoup
import time
from tqdm import tqdm, trange
import re
import pandas as pd




def crawler_modelseed(id):

    headers = {
        'user-agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.0.0 Safari/537.36",
        }
    url = 'https://modelseed.org/solr/compounds/select?wt=json&q=id:{}'.format(id)
    response = requests.request("GET", url, headers=headers)
    rawdata = BeautifulSoup(response.text, 'html.parser')
    try:
        abbreviation_r = re.findall(r'"abbreviation":".*"', str(rawdata))
        abbreviation = abbreviation_r[0].split('"')[3]
    except IndexError:
        abbreviation_r = abbreviation = 'nan'
    try:
        aliases_r = re.findall(r'"aliases":\["Name: .*', str(rawdata))
        aliases = aliases_r[0].split(';')[0].split(':')[2]
    except IndexError:
        aliases_r = aliases = 'nan'
    try:
        deltag_r = re.findall(r'"deltag":.*,', str(rawdata))
        deltag = deltag_r[0].split(':')[1].split(',')[0]
    except IndexError:
        deltag_r = deltag = 'nan'
    try:
        smiles_r = re.findall(r'"smiles":".*"', str(rawdata))
        smiles = smiles_r[0].split('"')[3]
    except IndexError:
        smiles_r = smiles = 'nan'
    try:
        kegg_r = re.findall(r'"KEGG: .*"', str(rawdata))
        kegg = kegg_r[0].split(' ')[1].split('"')[0].split(';')[0]
    except IndexError:
        kegg_r = kegg = 'nan'
    print('id:{} aliases:{}'.format(id, aliases))
    return [aliases, abbreviation, deltag, smiles, kegg, aliases_r, abbreviation_r, deltag_r, smiles_r, kegg_r, rawdata]



num = list(range(1, 15))
idlist = []
for n in num:
    if n <= 9:
        id = 'cpd' + '0000' + str(n)
    elif 10 <= n <= 99:
        id = 'cpd' + '000' + str(n)
    elif 100 <= n <= 999:
        id = 'cpd' + '00' + str(n)
    elif 1000 <= n <= 9999:
        id = 'cpd' + '0' + str(n)
    else:
        id = 'cpd' + str(n)
    idlist.append(id)

out = pd.DataFrame(index=idlist, columns=['aliases', 'abbreviation', 'deltag', 'smiles', 'kegg', 'aliases_r', 'abbreviation_r', 'deltag_r', 'smiles_r', 'kegg_r', 'rawdata'])
for i in trange(len(idlist)):
    time.sleep(0.25)
    data = crawler_modelseed(idlist[i])
    out.loc[idlist[i], :] = data

out2 = pd.DataFrame(index=idlist, columns=['aliases', 'abbreviation', 'deltag', 'smiles', 'kegg'])
out2.loc[:, 'aliases'] = out.loc[:, 'aliases']
out2.loc[:, 'abbreviation'] = out.loc[:, 'abbreviation']
out2.loc[:, 'deltag'] = out.loc[:, 'deltag']
out2.loc[:, 'smiles'] = out.loc[:, 'smiles']
out2.loc[:, 'kegg'] = out.loc[:, 'kegg']
out2.to_excel('./seeddata.xlsx')

