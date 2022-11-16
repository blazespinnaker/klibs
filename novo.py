import requests, os
import numpy as np, pandas as pd
import os, re
import matplotlib.pyplot as plt

def get_struct(seq, fname):
    if os.path.exists(fname):
        return
    esmfold_api_url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    r = requests.post(esmfold_api_url, data=seq)
    
    if r.status_code == 200:
        structure = r.text
        print (f'Prediction for is now complete and ready for download')
    else:
        print (r.status_code, 'Error for', seq,)
        structure = None
    with open(fname, 'w') as f:
        f.write(structure)


def get_dfa(only = []):
    d = None
    t = 0
    dfa = pd.DataFrame()
    for dirname, _, filenames in os.walk('/kaggle/input'):
        for filename in filenames:
            pth = os.path.join(dirname, filename)
            scre = re.compile('.*\/([0-9]{4}).csv', re.IGNORECASE).match(pth)
            if scre:
                sc = int(scre.group(1))
                if sc in only or True:
                    sc = sc / 10000
                    print(f"Adding {filename} to blend with score", sc)
                    df = pd.read_csv(pth)
                    dfa[scre.group(1)] = df['tm'].rank()
                    r = df['tm'].rank()
                    w = np.exp(11*sc)
                    r = r*w
                    t = t+w
                    d = r if d is None else d + r
                else:
                    print("Matched, but not in only",only, " - ", sc, pth)
            else:
                pass
                #print("Ignoring ", pth)
    dfa = dfa.sort_index(axis=1, ascending = False)
    with pd.option_context('display.max_rows', 2420, 'display.max_columns', None): 
        display(dfa.corr(method='spearman').style.background_gradient(cmap='coolwarm', axis=None).format(precision=3))
    return dfa
