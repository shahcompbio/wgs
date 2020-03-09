# Whole Genome Pipelines


Welcome to the home page for the whole genome sequencing pipelines documentation.


## Quick Setup

### Test Data
Test datasets for each of the five wgs subpipelines can be downloaded from azure storage. 
Just copy and paste the below commands for the subpipeline you which to run.
#### Alignment
```
wget https://wgstestsets.blob.core.windows.net/testsets/alignment.tar && tar -xvf alignment.tar 
rm alignment.tar && cd alignment
```
#### Realignment
```
wget https://wgstestsets.blob.core.windows.net/testsets/realignment.tar && tar -xvf realignment.tar 
rm realignment.tar && cd realignment
```
#### Variant Calling
```
wget https://wgstestsets.blob.core.windows.net/testsets/variant_calling.tar && tar -xvf variant_calling.tar 
rm variant_calling.tar && cd variant_calling
```
#### SV Calling
```
wget https://wgstestsets.blob.core.windows.net/testsets/sv_calling.tar && tar -xvf sv_calling.tar 
rm sv_calling.tar && cd sv_calling
```
#### Copy Number Calling
```
wget https://wgstestsets.blob.core.windows.net/testsets/cna.tar && tar -xvf cna.tar
rm cna.tar && cd cna
```
At this point, you should be in a new working directory with all the testdata you need to run the pipeline.
In addition, you will have three `yaml` files: `input.yaml`, `config.yaml` and `context.yaml`. 
<br/><br/>`input.yaml`: specifies the testdata files you are using.
<br/><br/>`config.yaml`: specifies subpipeline specific parameters.
<br/><br/>`context.yaml`: specifies information used by the pipeline to run locally in docker.

You will also have one or two python files: `make.py` and `make_ref.py`. If the tar you downloaded 
does not include `make_ref.py`, it just means all the reference data you need is included in your download.

### Customizing the yamls with local paths

`make.py` simply customizes the yaml files to use your working directory. Simply execute 
```
python make.py
```
If you do `make_ref.py` is included in the tar you downloaded, it means additional reference data is needed to 
run the subpipeline. This data was not included in the tar because of its large size. Go to __TODO: GET REF DATA PAGE___
To download the reference data to run the given subpipeline. Then, execute
```
python make_ref.py {directory with downloaded reference data}
```
Once you've run these two python files, your ready to make the pipeline run command.

### Running a Subpipeline with Testdata locally

To run the pipeline locally, you'll need to construct two files that combine the `yaml` inputs with subpipeline-specific commands.
The first file is called `final_run.sh`, and it defines the docker image you are using to run the pipeline along with docker-speficic parameters. The second file, `run_pipeline.sh`, contains the command that will run the pipeline and is passed into `final_run.sh`. 

#### `final_run.sh`:
Create a file named `final_run.sh` and copy the following into it.
By default, use the `wgspipeline/wgs:v0.0.4` docker image to run a subpipeline. 
```
docker run -v $PWD:$PWD -w $PWD -v /var/run/docker.sock:/var/run/docker.sock \
  -v `which docker`:`which docker` wgspipeline/wgs:v0.0.4 sh run_pipeline.sh
```
You'll notice that this command uses `run_pipeline.sh`.

Note: If you want to test out the docker image before running the pipeline, simply run
```
docker run -it -v $PWD:$PWD -w $PWD -v /var/run/docker.sock:/var/run/docker.sock \
  -v `which docker`:`which docker` wgspipeline/wgs:v0.0.4 bash
which wgs
```
#### `run_pipeline.sh`:

`run_pipeline.sh` contains the actual command that runs the wgs pipeline, so it is different for every subpipeline. 

For the given subpipeline you want to run, copy and paste below into `run_pipeline.sh`.

#### Alignment
```
  wgs alignment --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local --config_file config.yaml \
  --context_config context.yaml
```
#### Realignment
```
  wgs realignment --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local --config_file config.yaml \
  --context_config context.yaml
```
#### Variant Calling
```
  wgs variant_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local --config_file config.yaml \
  --context_config context.yaml
```
#### Breakpoint Calling
```
  wgs breakpoint_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local --config_file config.yaml \
  --context_config context.yaml
```
#### Copy Number Calling
```
  wgs copynumber_calling --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local --config_file config.yaml \
  --context_config context.yaml
```

Once both `final_run.sh` and `run_pipeline.sh` are made, simply run 
```
sh final_run.sh
```
To launch the pipeline.
