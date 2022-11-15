import requests, os
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
