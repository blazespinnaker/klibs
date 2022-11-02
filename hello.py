from kaggle_secrets import UserSecretsClient
user_secrets = UserSecretsClient()
gcloud_api = user_secrets.get_secret("gcloud")
print(gcloud_api, file=open("/tmp/key.json", "w"))
!gcloud auth activate-service-account --key-file /tmp/key.json
!mkdir data
!gsutil -m cp gs://nesp/data.zip data
!cd data && unzip -q data.zip
!rm data/data.zip
