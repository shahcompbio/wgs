
  
  
  
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
tar -xvf alignment_data.tar.gz
cd data
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
  --refdir ref --maxjobs 4
```

Please refer to the docker guide to learn how to launch pipeline.

#### Variant Calling
1. download test datasets
```
wget https://wgstestsets.blob.core.windows.net/datasets/variant_data.tar.gz
tar -xvf variant_data.tar.gz
cd data
```
2. create input.yaml
```
HCC1395:
  normal: data/HCC1395BL_chr15_snps.bam
  normal_id: HCC1395BL
  tumour: data/HCC1395_chr15_snps.bam
  tumour_id: HCC1395
```
3. create pipeline.sh
```
#!/bin/bash

wgs variant_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir ref --maxjobs 4
```

### Copynumber calling
1. download test datasets
```
wget https://wgstestsets.blob.core.windows.net/datasets/copynumber_data.tar.gz
tar -xvf copynumber_data.tar.gz
cd data
```
2. create input.yaml
```
HCC1395:
  normal: data/HCC1395BL_chr15_snps.bam
  normal_id: HCC1395BL
  tumour: data/HCC1395_chr15_snps.bam
  tumour_id: HCC1395
  target_list: data/targets.tsv
```
3. create pipeline.sh
```
#!/bin/bash

wgs copynumber_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir ref --titan --maxjobs 4
```

### Breakpoint calling
1. download test datasets
```
wget https://wgstestsets.blob.core.windows.net/datasets/breakpoint_data.tar.gz
tar -xvf breakpoint_data.tar.gz
cd data
```
2. create input.yaml
```
SA123:
  normal: data/normal.bam
  normal_id: SA123N
  tumour: data/small.bam
  tumour_id: SA123
```
3. create pipeline.sh
```
#!/bin/bash

wgs breakpoint_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir ref --maxjobs 4
```

###  Realignment
1. download test datasets
```
wget https://wgstestsets.blob.core.windows.net/datasets/realignment_data.tar.gz
tar -xvf realignment_data.tar.gz
cd data
```
2. create input.yaml
```
SA123:
  input: data/A20875_3_lanes_dupsFlagged_chr22_paired.bam
  output: bams/output.bam
```
3. create pipeline.sh
```
#!/bin/bash

wgs realignment --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir ref --maxjobs 4
```

###  Postprocessing
1. download test datasets
```
wget https://wgstestsets.blob.core.windows.net/datasets/postprocessing_data.tar.gz
tar -xvf postprocessing_data.tar.gz
cd data
```
2. create input.yaml
```
Sample_123:
  normal: data/bams/normal.bam
  tumour: data/bams/variants.bam
  variant_dir: data/variants
  breakpoint_dir: data/breakpoints
  copynumber_dir: data/copynumber
```
3. create pipeline.sh
```
#!/bin/bash

wgs postprocessing --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local \
  --refdir ref --maxjobs 4
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

docker run --rm -v $PWD:$PWD -w $PWD  -v /var/run/docker.sock:/var/run/docker.sock -v `which docker`:`which docker` wgspipeline/wgs:v0.0.6 sh pipeline.sh
```

4. start the pipeline:
```
sh launcher.sh
```