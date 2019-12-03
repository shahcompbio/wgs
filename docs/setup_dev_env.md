
#### Requirements:

this guide assumes that you have
1. write access to wgspipeline dockerhub org
2. docker installed on the node
3. root access


# Setup


### create working dir and clone code

```
mkdir /devsrc
cd /devsrc
git clone https://github.com/shahcompbio/wgs.git
```

### build docker container

save the following in `dockerfile` file in `/devsrc`
```
# build on top of out base image
FROM wgspipeline/python_wgs:v0.0.1

# Install any needed packages specified in requirements.txt
RUN rm -rf /opt/conda/lib/python2.7/site-packages/pypeliner* /opt/conda/lib/python2.7/site-packages/wgs* /opt/conda/lib/python2.7/site-packages/biowrappers*
RUN pip install git+https://github.com/shahcompbio/pypeliner.git@master
RUN pip install git+https://bitbucket.org/aroth85/biowrappers.git@singlecell
RUN pip install dill
# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME wgs

ENV PYTHONPATH /devsrc/wgs

# Run app.py when the container launches
CMD ["wgs"]
```

build:
`docker build -t wgspipeline/wgs:v0.0.1 .`

push:
`docker push wgspipeline/wgs:v0.0.1`


### Context config

save the following yaml to a file.


```
docker:
    server: 'docker.io'
    username: null
    password: null
    mounts:
      refdata: /refdata
      datadrive: /devsrc
```


##### reference data
You need to download the reference data to `/refdata' directory on the node. The data can be downloaded from blob. 

the data is in `wgscomputedev` storage account in `refdata` container.


#### Alignment pipeline:

1. test data:

The input yaml should have the following format.
```
SA123:
  fastqs:
    L001:
      fastq1: /devsrc/testdata/R1.fastq.gz
      fastq2: /devsrc/testdata/R2.fastq.gz
  bam: results/bams/SA123.bam
```


`R1.fastq.gz` and `R2.fastq.gz` can be found in blob under wgscomputedev storage account and alignmenttestdata container. Please contact <grewald@mskcc.org> for access.


save the following in a shell script:

```
cd /devsrc/wgs && python setup.py develop && cd /devsrc

wgs alignment \
--input_yaml input.yaml \
--out_dir output \
--tmpdir temp \
--pipelinedir pipeline \
--loglevel DEBUG \
--submit local \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' \
--sentinel_only  --maxjobs 100 \
--config_override '{"cluster":"juno"}' \
--context_config path_to_context_config_yaml_here
```

and then run the script:

```

docker run -w $PWD -v $PWD:$PWD --rm -v /var/run/docker.sock:/var/run/docker.sock \
-v /usr/bin/docker:/usr/bin/docker  wgspipeline/wgs:v0.0.2 bash run.sh
```


