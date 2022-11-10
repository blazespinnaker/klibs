from kaggle_secrets import UserSecretsClient
user_secrets = UserSecretsClient()
gcloud_api = user_secrets.get_secret("gcloud")
print(gcloud_api, file=open("/tmp/key.json", "w"))
!gcloud auth activate-service-account --key-file /tmp/key.json
!mkdir data
!gsutil -m cp gs://nesp/data.zip data
!cd data && unzip -q data.zip
!rm data/data.zip

!zip pdbs.zip  -r ./

import os
#cmd = f"kaggle datasets version -p /kaggle/working/optuna -m 'new_{count}'"
cmd = f"kaggle datasets init -p /kaggle/working"
os.system(f"export KAGGLE_USERNAME=kaggleqrdl;export KAGGLE_KEY=225a4247dea6373cf03067ba0965d995;{cmd}")

tm5210 = dfa['6010'].copy()
for i,r in novo.iterrows():
    tm5210.loc[r['novo_ind']] = r['dtm_pred']
print(len(tm5210))
dfa.corrwith(tm5210,method = 'spearman')

samp = pd.read_csv("../input/novozymes-enzyme-stability-prediction/sample_submission.csv")
samp['tm'] = pd.Series(tm5210)#*0.05 + dfa['6010'].rank()
samp.set_index("seq_id").to_csv("pos_cor.csv")
display(dfa.corrwith(samp['tm'], method= 'spearman'))



from sklearn.model_selection import cross_val_score
def fs(clff, df, Xcols, ycol, basep = []):
    base = basep.copy()
    metamaxv = -999
    for i in range(len(Xcols)):
        maxv = -999
        for c in Xcols:
            if c in base:
                continue
            base_new = base + [c]
            model = clff(0)
            #score = cross_val_score(model,df[base_new], df[ycol] )
            model.fit(df[base_new], df[ycol] )
            prds = model.predict(novo.loc[:,base_new])
            novo['dtm_pred'] = prds
            retval = novo['dtm'].corr(novo['dtm_pred'], method = "spearman")
            #retval = np.mean(score)
            if retval > maxv:
                maxv = retval
                cols = base_new
                print(maxv, cols)
        base = cols
        if maxv > metamaxv:
            metamaxv = maxv
            metacols = cols
            print("new meta", maxv, cols)
    return metacols, metamaxv
