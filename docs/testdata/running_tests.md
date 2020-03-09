# WGS test data

## Realignment

docker command line:

```
docker run -v $PWD:$PWD -w $PWD -v /var/run/docker.sock:/var/run/docker.sock \
  -v `which docker`:`which docker` wgspipeline/wgs:v0.0.4 \
  wgs realignment --input_yaml input.yaml \
  --out_dir output --tmpdir temp --pipelinedir pipeline \
  --loglevel DEBUG --submit local --config_file config.yaml \
  --context_config context.yaml
```

## Alignment

```
docker run -v $PWD:$PWD -w $PWD -v /var/run/docker.sock:/var/run/docker.sock \
  -v `which docker`:`which docker` wgspipeline/wgs:v0.0.4 \
  wgs alignment --input_yaml alignment/input.yaml \
  --out_dir alignment/output --tmpdir alignment/temp --pipelinedir alignment/pipeline \
  --context_config context.yaml --config_file alignment/config.yaml \
  --loglevel DEBUG --submit local
```

## Variant Calling

```
docker run -v $PWD:$PWD -w $PWD -v /var/run/docker.sock:/var/run/docker.sock \
  -v `which docker`:`which docker` wgspipeline/wgs:v0.0.4 \
  wgs variant_callinlsg --input_yaml variant_calling/input.yaml \
  --out_dir variant_calling/output --tmpdir variant_calling/temp --pipelinedir variant_calling/pipeline \
  --context_config context.yaml \
  --loglevel DEBUG --submit local
```
 
 ## Breakpoint Calling

```
docker run -v $PWD:$PWD -w $PWD -v /var/run/docker.sock:/var/run/docker.sock \
  -v `which docker`:`which docker` wgspipeline/wgs:v0.0.4 \
  wgs breakpoint_calling --input_yaml sv_calling/input.yaml \
  --out_dir sv_calling/output --tmpdir sv_calling/temp --pipelinedir sv_calling/pipeline \
  --context_config context.yaml \
  --loglevel DEBUG --submit local
```
 