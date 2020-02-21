# wgs on a  cluster 

#### Requirements:

1. singularity install on head node and all compute nodes.

      Use ```module load singularity``` on `juno` to activate singularity.
      
2. System should allow ssh into localhost without any passphrase or other issues on head node and all compute nodes. The following command should work on all nodes
   ```
   ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null localhost
   ```
    this can be done through ssh keys on the cluster. We don't support ssh keys with a passphrase at the moment.
3. access to dockerhub on all nodes
    ```
    singularity run docker://docker.io/singlecellpipeline/single_cell_pipeline:v{VERSION}
    ```

### Download the reference data 

The pipeline reference data is available at the following locations:

*Juno cluster at MSKCC:*

```
/juno/work/shah/reference/wgs_pipeline
```

*Shahlab cluster at GSC:*
```
/shahlab/pipelines/reference/wgs_pipeline
```

*Azure:*
```
https://wgscomputedata.blob.core.windows.net/referrence-*
```

NOTE: Access to data on azure is restricted to the shahlab. To request access to the data please contact Diljot Grewal <grewald@mskcc.org.>

### Download the test data

see `/docs/testdata` for instructions on downloading testdata and test input yamls.

### Create the context configuration file

The context configuration file contains details about the docker container registry and the directories that need to be mounted inside the singularity containers. Optionally you can also define some default job parameters. The name_match can be used to selectively apply the parameters to the jobs based on their name.

The following example points the pipeline to the dockerhub container registry. The mounts section lists the directories that must be accessible to the pipeline. This includes the directories that contain inputs, output, reference data and input files such as the yaml inputs.
the context section of the yaml snippet below sets the walltime and number of times the job will be retried. The walltime parameter will be set when we launch the pipeline later. 

```
singularity:
    server: 'docker.io'
    username: null
    password: null
    local_cache: '/juno/work/shah/runner/singularity/cache'
    singularity_exe: 'singularity'
    mounts:
      juno: /juno/work/shah/pipelinedir
      reference: /juno/work/shah/reference
      common: /common
context:
  alljobs:
    name_match: '*'
    ctx:
      walltime: '4:00'
      walltime_num_retry: 5
      walltime_retry_increment: '48:00'

```
Please update the local_cache and the mounts directories before you run. The `\common` mount will mount the required LSF paths in singularity to give job submission access to the pipeline. If the LSF executables are located at a different location, please provide that path in the config.


### launch the pipeline

write the following to a file:

```
export PATH=/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin:$PATH

wgs alignment --input_yaml input.yaml --out_dir output --tmpdir temp \
--pipelinedir pipeline --loglevel DEBUG --submit lsf \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' 
--sentinel_only  --maxjobs 100 --config_override '{"cluster":"juno"}' 
--context_config context_config.yaml --nocleanup
```

Please refer to [doc](../../README.md) for detailed instructions for running all single cell sub commands.

launch the pipeline:

```
export SINGULARITY_CACHEDIR=/juno/work/shah/runner/singularity/cache
singularity run --bind /common --bind /juno/work  docker://docker.io/singlecellpipeline/wgs_pipeline:v{VERSION} sh /path/to/shell/script/from/previous/step
```

The `--bind /common` will mount the `/common` directory inside the singularity. The PATH environment variable must also be set to point to the location of LSF binaries. This will make the commands such as `bsub`, `bjobs` and `bhosts` available to the pipeline. This path will depend on the singularity location. 


This command runs the alignment portion of the wgs pipeline. See [doc](../../README.md) for details on inputs to other wgs pipelines.
