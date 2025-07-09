This doc is for setting up Google batch processing on the all of us platform. In the end a local execution of the pipeline was performed because staging the VEP cache to the VM was a bottleneck. 

### Set-Up VEP Container
1. Download docker image on local:
```bash
docker pull ensemblorg/ensembl-vep:release_113.4

```
2. Tag and upload to private Google artifact registry:
```bash
docker tag ensemblorg/ensembl-vep:release_113.4 us-central1-docker.pkg.dev/cnv-profiler-docker/allofus-jacq/ensemblorg/ensembl-vep:release_113.4

docker push us-central1-docker.pkg.dev/cnv-profiler-docker/allofus-jacq/ensemblorg/ensembl-vep:release_113.4
```

### Configure Google Credentials

1. On the VM, we need to validate our credentials and update our google configurations:

```bash
gcloud auth application-default login
```
2. Next we need to update a nextflow variable to point to the google configuration file:

```bash
export GOOGLE_APPLICATION_CREDENTIALS='/home/jupyter/.config/gcloud/application_default_credentials.json'
```
3. Nextflow can complain if the project ID is not found in the credentials file. Make sure it looks like this:

```json
{
  "account": "",
  "project_id": "terra-vpc-sc-2393c802",
  "client_id": "<redacted>",
  "client_secret": "<redacted>",
  "quota_project_id": "terra-vpc-sc-2393c802",
  "refresh_token": "<redacted>",
  "type": "authorized_user",
  "universe_domain": "googleapis.com"
}
```

### Update Config

Nextflow will require a specific configuration to enable batch processing. Follow instructions here [https://www.nextflow.io/docs/latest/google.html](https://www.nextflow.io/docs/latest/google.html): 

The process using google batch should look like this:

```
docker.enabled = true

process VEP {
	executor = 'google-batch'
	container = 'docker://us-central1-docker.pkg.dev/cnv-profiler-docker/allofus-jacq/ensemblorg/ensembl-vep'
	location = 'us-central1'
}

```

If starting resources are found in a bucket,  a -bucket-dir parameter can be specified in the command line options. 


