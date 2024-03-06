'''
@qqlaoxia 20240306
'''

import argparse
import requests
import pandas as pd
import json
import os
import time
def extract_data(i):
    cmp='YMDB'+'{:05d}'.format(i)
# Send an HTTP request to get the content of the web page
    url = r'https://www.ymdb.ca/compounds/'+cmp
    response = requests.get(url)
    web_content = response.text
    data = json.loads(web_content)

    if isinstance(data, dict):
        data = [data]
    extracted_data = []
    for item in data:
        row_data = {key: value for key, value in item.items()}
        extracted_data.append(row_data)

    df = pd.DataFrame(extracted_data)
    # if not os.path.exists(cmp):
    #     os.mkdir(cmp)
    #
    # # Save the DataFrame to an Excel file
    # excel_path = './%s/%s.xlsx' %(cmp,cmp)
    # df.to_excel(excel_path, index=False)
    # print(f'Data saved to {excel_path}')


    # % rearrange data
    df_transposed = df.T.reset_index()
    df_transposed.columns = ['Key', 'Content']

    # Save the transposed DataFrame to a new Excel file
    transposed_excel_path = f'../data/YMDB/{cmp}_transposed.csv'
    df_transposed.to_csv(transposed_excel_path, index=False, header=True)
    print(f'Transposed data saved to {transposed_excel_path}')
    time.sleep(5)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='input number',type=int)
    args = parser.parse_args()
    extract_data(args.i)
    # It is suggested to used parallel arithmetic in matlab.
    # parfor i=1:16338
    # system(['python get_YMDB_data.py -i ',num2str(i)]);
    # end



