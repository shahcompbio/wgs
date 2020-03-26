
# Whole Genome Pipelines


Welcome to the home page for the whole genome sequencing pipelines documentation.


## Quick Setup

### Test Data
Test datasets for each of the five wgs subpipelines can be downloaded from azure storage. 
Just copy and paste the below commands for the subpipeline you which to run.

#### Alignment

1. Download the test data set:
```
wget  https://wgstestsets.blob.core.windows.net/datasets/alignment_data.tar.gz
```
2. create input.yaml
```
TEST:
  fastqs:
    simulated_lane:
      fastq1: data/test.r1.fastq
      fastq2: data/test.r2.fastq
  bam: bams/TEST.bam
  readgroup_info:
    ID: '{lane_id}'
    PU: '{lane_id}'
    SM: '{sample_id}'
    LB: 'TEST'
    CN: 'BCCAGSC'
    PL: 'ILLUMINA'
```
3. create pipeline.sh file: 
```
#!/bin/bash

wgs alignment --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir ref
```

Please refer to the docker guide to learn how to launch pipeline.

#### Variant Calling
1. download test datasets
```
wget https://wgstestsets.blob.core.windows.net/datasets/variant_data.tar.gz
tar -xvf variant_data.tar.gz
```
2. create input.yaml
```
HCC1395:
  normal: data/data/HCC1395BL_chr15_snps.bam
  normal_id: HCC1395BL
  tumour: data/data/HCC1395_chr15_snps.bam
  tumour_id: HCC1395
```
3. create pipeline.sh
```
#!/bin/bash

wgs variant_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir data/ref
```

#### Run with docker

1. create context.yaml file:
```
docker:
    server: 'docker.io'
    org: wgspipeline
    username: null
    password: null
    mounts:
      local: <your current working dir>
```

2. add context.yaml file to pipeline.sh. For instance, for variant calling
```
wgs variant_calling --input_yaml inputs.yaml
...
--context_config context.yaml
```

3. create launcher.sh
```
#!/bin/bash

docker run -it -v $PWD:$PWD -w $PWD  -v /var/run/docker.sock:/var/run/docker.sock -v `which docker`:`which docker` wgspipeline/wgs:v0.0.6 sh pipeline.sh
```

4. start the pipeline:
```
sh launcher.sh
```
