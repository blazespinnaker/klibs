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


def get_dfa(only = [], dirf = "/home/tim/subs"):
    d = None
    t = 0
    dfa = pd.DataFrame()
    for dirname, _, filenames in os.walk(dirf):
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


import pandas as pd, numpy as np, nltk
#size of test is 2413
def avg_sub(sub, use_test = None):
    wt = "VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK"
    if use_test is None:
        test_df = pd.read_csv("/kaggle/input/novozymes-enzyme-stability-prediction/test.csv")
        test_df.loc[1169, 'protein_sequence'] = wt[:-1]+"L" #1169 is the same as wt
        test_df['wt_sequence'] = wt
    else:
        test_df = use_test
    test_df['edit_idx'] = test_df.apply(lambda x:[i for i in range(len(x['protein_sequence'])) if x['protein_sequence'][i] != x['wt_sequence'][i]][0], axis = 1)
    test_df['mttype'] = test_df.apply(lambda x:"deletion" if len(x['protein_sequence']) == 220 else "substitution", axis = 1)
    test_df['mt'] = test_df.apply(lambda x:x['protein_sequence'][x['edit_idx']], axis = 1)
    test_df['wt'] = test_df.apply(lambda x:x['wt_sequence'][x['edit_idx']], axis = 1)
    test_df['wtmt'] = test_df.apply(lambda x:x['wt']+str(x['edit_idx']+1)+x['mt'] if x['mttype'] == "substitution" else x['wt']+str(x['edit_idx']+1)+"_", axis = 1)
    test_df['fasta'] = test_df.apply(lambda x:"fungi_"+x['wt']+str(x['edit_idx']+1)+x['mt'], axis = 1)
    test_df['maestro'] = test_df.apply(lambda x:x['wt']+str(x['edit_idx']+1)+".A{"+x['mt']+"}", axis = 1)
    test_df['fasta'] = test_df.apply(lambda x:"fungi_"+x['wt']+str(x['edit_idx']+1)+x['mt'], axis = 1)
    test_df['tm'] = sub.rank()
    means = test_df.query("mttype == 'substitution'").sort_values('edit_idx').groupby('edit_idx').mean()['tm']
    try:
        test_df['edit_idx_avg'] = test_df.apply(lambda x:means.loc[x['edit_idx']], axis = 1)
    except:
        return means
    test_df['edit_idx_avg_plus_minus'] = test_df.apply(lambda x:(x['tm']-x['edit_idx_avg'])/x['edit_idx_avg'], axis = 1)
    test_df['edit_idx_avg_plus_minus_const'] = test_df.apply(lambda x:(x['tm']-x['edit_idx_avg']), axis = 1)
    #test_df['edit_idx_avg_plus_minus_avg'] = test_df.apply(lambda x:x['edit_idx_avg']-x['edit_idx_avg'], axis = 1)
    mt_means = test_df.query("mttype == 'substitution'").groupby('mt').mean()['edit_idx_avg_plus_minus']
    mt_means_const = test_df.query("mttype == 'substitution'").groupby('mt').mean()['edit_idx_avg_plus_minus_const']
    try:
        test_df['edit_idx_avg_plus_minus_mt'] = test_df.apply(lambda x:mt_means.loc[x['mt']] , axis = 1)
    except:
        return mt_means
    test_df['edit_idx_avg_plus_minus_adj'] = test_df.apply(lambda x:(x['edit_idx_avg_plus_minus']*4 + x['edit_idx_avg_plus_minus_mt']*1)/5 , axis = 1) #try weighting these two, eg x['edit_idx_avg_plus_minus']*0.75
    test_df['edit_idx_avg_plus_minus_mt_const'] = test_df.apply(lambda x:mt_means_const.loc[x['mt']] , axis = 1)
    test_df['edit_idx_avg_plus_minus_adj_const'] = test_df.apply(lambda x:(x['edit_idx_avg_plus_minus_mt_const']*0.5) , axis = 1) #try weighting these two, eg x['edit_idx_avg_plus_minus']*0.75
    test_df['adj_tm'] =  test_df.apply(lambda x:x['edit_idx_avg']*(1+x['edit_idx_avg_plus_minus_adj']) , axis = 1)
    test_df['adj_tm_new_rank'] = test_df['adj_tm'].rank()
    test_df['adj_tm_const'] =  test_df.apply(lambda x:x['tm']+x['edit_idx_avg_plus_minus_adj_const'] , axis = 1)
    test_df['adj_tm_new_rank_const'] = test_df['adj_tm_const'].rank()
    #display(test_df.corr())
    test_df['orig_tm'] = test_df['tm']
    return test_df

def fs(params):
    use_base = params['use_base']
    selcb = params['selcb']
    enc = params['enc']
    maxcb = params['maxcb']
    metacb = params['metacb']
    start_base = params['start_base']
    nnovoc = params['nnovoc']
    if params['ycol'] in enc:
        enc = [x for x in enc if x != params['ycol']]
    if start_base != []:
        print("Start base",selcb(start_base, enc, params, nnovoc))
    metamaxv = -999
    reversemax = -999
    cols = []
    metacols = []
    base = start_base
    for i in range(len(use_base)+1):
        maxv = -999
        for p in use_base: 
            if p in base:
                continue
            base_new = base + [p]
            retval, model = selcb(base_new, enc, params, nnovoc)
            if retval > maxv:
                forp = p
                maxv = retval
                cols = base_new
                maxcb("forward", retval,cols, model,enc,params)
        base = cols
        if len(base) > 1:
            for p in base:
                if p == forp:
                    continue
                base_new = [x for x in base if x != p]
                retval, info = selcb(base_new, enc, params, nnovoc)
                if retval > reversemax and retval > maxv:
                    reversemax = retval
                    cols = base_new
                    maxcb("\nreverse", retval, cols, model,enc,params)
        base = cols
        if maxv > metamaxv:
            metamaxv = maxv
            metacols = cols
            metacb("new", i, metamaxv, metacols)
        else:
            metacb("still", i, metamaxv, metacols)

