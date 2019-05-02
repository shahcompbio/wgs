# Whole Genome Pipelines




### setup and Installation:

Set up conda with the required packages.



### From Source
Add channels:
```
conda config --add channels shahcompbio
conda config --add channels dranew
conda config --add channels aroth85
conda config --add channels componc
conda config --add channels bioconda
```

Then create an environment with the required packages:

```
conda create --name wgspipeline --file conda_packages.txt
```

Activate the environment:

```
source activate wgspipeline
```

Install pacakges from source:

```
pip install git+https://bitbucket.org/aroth85/biowrappers.git@singlecell
pip install git+https://github.com/shahcompbio/pypeliner.git@master
pip install git+https://github.com/shahcompbio/wgs.git@master
pip install git+https://dgrewal@svn.bcgsc.ca/bitbucket/scm/~dgrewal/vizutils.git
pip install git+https://dgrewal@svn.bcgsc.ca/bitbucket/scm/museq/museqportrait.git@new_vcf_format
```

#### Input File format

```
SAMPLE_ID:
  fastqs:
    normal:
      NORMAL_SAMPLE_LANE_1_ID:
        fastq1: /path/to/fastq_r1.fastq.gz
        fastq2: /path/to/fastq_r2.fastq.gz
      NORMAL_SAMPLE_LANE_2_ID:
        fastq1: /path/to/fastq_r1.fastq.gz
        fastq2: /path/to/fastq_r2.fastq.gz
    tumour:
      TUMOUR_SAMPLE_LANE_1_ID:
        fastq1: /path/to/fastq_r1.fastq.gz
        fastq2: /path/to/fastq_r2.fastq.gz
      TUMOUR_SAMPLE_LANE_2_ID:
        fastq1: /path/to/fastq_r1.fastq.gz
        fastq2: /path/to/fastq_r2.fastq.gz
      TUMOUR_SAMPLE_LANE_3_ID:
        fastq1: /path/to/fastq_r1.fastq.gz
        fastq2: /path/to/fastq_r2.fastq.gz
  normal: /path/to/output/aligned/normal.bam
  normal_id: NORMAL_SAMPLE_ID
  tumour: /path/to/output/aligned/tumour.bam
  tumour_id: TUMOUR_SAMPLE_ID
  breakpoints: /path/to/destruct/breakpoints.csv
```

The fastqs section is only required for the alignment workflow and the full workflow (if the alignment flag is set).
The breakpoints section is only required for the copynumber workflow if you need remixt results. 

#### Launch Full workflow

```
wgs all --input_yaml input.yaml --out_dir results --tmpdir tmp --pipelinedir pipeline --submit lsf --maxjobs 1000 --nocleanup --loglevel DEBUG --nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'  --config_override '{"cluster":"juno"}' --context_config context.yaml --alignment --sentinal_only --rerun
```

#### Launch Alignment workflow

```
wgs alignment --input_yaml input.yaml --out_dir results --tmpdir tmp --pipelinedir pipeline --submit lsf --maxjobs 1000 --nocleanup --loglevel DEBUG --nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'  --config_override '{"cluster":"juno"}' --context_config context.yaml --alignment --sentinal_only --rerun
```

#### Launch variant calling workflow

```
wgs variant_calling --input_yaml input.yaml --out_dir results --tmpdir tmp --pipelinedir pipeline --submit lsf --maxjobs 1000 --nocleanup --loglevel DEBUG --nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'  --config_override '{"cluster":"juno"}' --context_config context.yaml --alignment --sentinal_only --rerun
```



#### Launch copynumber calling workflow

```
wgs copynumber_calling --input_yaml input.yaml --out_dir results --tmpdir tmp --pipelinedir pipeline --submit lsf --maxjobs 1000 --nocleanup --loglevel DEBUG --nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'  --config_override '{"cluster":"juno"}' --context_config context.yaml --alignment --sentinal_only --rerun
```

#### Launch breakpoint calling workflow

```
wgs breakpoint_calling --input_yaml input.yaml --out_dir results --tmpdir tmp --pipelinedir pipeline --submit lsf --maxjobs 1000 --nocleanup --loglevel DEBUG --nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'  --config_override '{"cluster":"juno"}' --context_config context.yaml --alignment --sentinal_only --rerun
```


#### Common Options:


##### Submit

`--submit lsf` to run on LSF clusters
`--submit local` to run locally
`--submit asyncqsub` to run on SGE based cluster

##### nativespec

use `--nativespec` to specify the cluster job submission format. You can use the following keywords as place holder and pipeline will automatically decide the best values for the jobs.

we support the following:
`{mem}` will be replaced with the optimal memory usage for each job
`{ncpus}` will be replaced with the optimal number of cpus for each job
`{walltime}` will be replaced with the optimal walltime for each job

These parameters will be passed to the job scheduler when running the pipeline.

For instance on a LSF based cluster, the nativespec might look like the following:

`--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"'`

##### sentinel only:
The pipeline looks at the files in the filesystem on reruns to track completed jobs. On some filesystems this might cause slowdowns. To replace this with a database please specify `--sentinel_only`


##### config options
The pipeline defaults to preset values for most configuration parameters. You can change these parameters by:

###### custom config file:
1. generate a new config file with

```
wgs generate_config --pipeline_config config.yaml
```
2. open the generated config yaml file, make changes where necessary and save it.
3. launch the pipeline with the `--config_file /path/to/config` parameter.

###### config override
You can also override certain values in the config file with the `--config_override` parameter. The config_override and config_file are mutually exclusive options.

The config override option accepts a json object. this json will override values in the internal config file. please generate a new config file for reference.

The pipeline also comes with some presets for config override. For instance:
1. if you're running this pipeline on MSKCC's juno cluster, please specify `--config_override '{"cluster":"juno"}'`.
2. If you're running the pipeline on BCCRC's shahlab cluster, please specify `--config_override '{"cluster":"shahlab"}'`

##### rerun
`--rerun` will run all jobs again, even if they've been run before.

##### context config

you can also specify a context config file to override the job execution parameters for certain job types. 
For instance:

`--context_config context.yaml` 
where context.yaml is
```
context:
  alljobs:
    name_match: '*'
    ctx:
      walltime: '04:00'
      walltime_num_retry: 5
      walltime_retry_increment: '48:00'
```
will update all jobs to 4 hrs of walltime and the pipeline will retry each job up to 5 times on failure and increment walltime by 2 days on each retry.

##### maxjobs
specifies the maximum number of jobs that pipeline will run in parallel.

##### nocleanup
do not clean up intermediates

##### loglevel
logging level. 

