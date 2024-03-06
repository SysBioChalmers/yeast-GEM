####################################################################
# get concentration
# The concentration range is calculated according to the following logic
# For example: In YMDB, three conc is recorded for glutamine ('17140.0 ± 1971.0 umol/L', '2500.0 ± 500.0 umol/L', '15000.0 ± 0.0 umol/L')
# maximal conc: The maximum concentration plus the corresponding error. In this case, is "17140.0 + 1971.0 umol/L"
# minimal conc: The minimum concentration minus the corresponding error. In this case, is "2500.0 - 500.0 umol/L"


import os
import numpy as np
import pandas as pd
import re


def get_conc(rawconc):
    concSplit = rawconc.split('}, ')
    # unit: umol/L
    all_conc = []
    # unit: mol/L
    maxmin = []
    for c in concSplit:
        regex = r'\'([^\']*)\''
        result = re.findall(regex, c)
        try:
            if (result[5] == '&#181;M') or (result[5] == 'uM') or (result[5] == 'umol/L'):
                all_conc.append(str(result[3] + ' ± ' + result[7] + ' umol/L'))
                maxmin.append((float(result[3]) - float(result[7]))*10**(-6))
                maxmin.append((float(result[3]) + float(result[7]))*10**(-6))
        except IndexError:
            pass
    return all_conc, maxmin

# load data and put them together
folder_path = r"../data/YMDB"
fileName = os.listdir(folder_path)
data = pd.DataFrame()
for fn in fileName:
    file_path = os.path.join(folder_path, fn)
    file = pd.read_excel(file_path)
    id = file.iloc[0, 1]
    if id == 'YMDB00001':
        data.index = list(file.loc[:, 'Key'])[1:]
    else:
        pass
    data.loc[:, id] = list(file.loc[:, 'Content'])[1:]
    print(fn)
# data.to_csv(r'../data/allConcData.csv')

# compute maximal conc and minimal conc
conc = pd.DataFrame(index=['conc', 'maxConc', 'minConc'], columns=data.columns)
for i in range(len(data.columns)):
    print(i)
    rawconc = data.iloc[44, i]
    all_conc, maxmin = get_conc(rawconc)
    if len(all_conc) != 0:
        conc.loc['conc', data.columns[i]] = all_conc
        conc.loc['maxConc', data.columns[i]] = np.max(maxmin)
        conc.loc['minConc', data.columns[i]] = np.min(maxmin)
    else:
        pass

data = pd.concat([data, conc], axis=0)
newdata = data.loc[['name', 'chebi_id', 'kegg_id', 'conc', 'maxConc', 'minConc']]
cleaned_data = newdata.dropna(axis=1, subset=['conc'])
cleaned_data.to_excel('../data/allConcData.xlsx')
# # test unit
#
# def check_unit(rawconc):
#     concSplit = rawconc.split('}, ')
#     unit = []
#     try:
#         for c in concSplit:
#             regex = r'\'([^\']*)\''
#             result = re.findall(regex, c)
#             unit.append(result[5])
#     except:
#         pass
#     return unit
#
# all_unit= []
# for i in range(len(data.columns)):
#     rawconc = data.iloc[44, i]
#     unit = check_unit(rawconc)
#     for u in unit:
#         if u not in all_unit:
#             all_unit.append(u)
#             print(data.columns[i])
#             print(u)