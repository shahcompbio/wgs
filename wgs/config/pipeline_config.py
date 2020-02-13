import collections
import os
import warnings

import wgs
import yaml
from wgs.utils import helpers

from azure_config import azure_config
from juno_config import juno_config
from shahlab_config import shahlab_config


def get_version():
    version = wgs.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


def containers():
    version = get_version()
    docker_images = {
        'bwa': 'wgspipeline/bwa:v0.0.1',
        'samtools': 'wgspipeline/samtools:v0.0.1',
        'picard': 'wgspipeline/picard:v0.0.1',
        'wgs': 'wgspipeline/wgs:v{}'.format(version),
        'strelka': 'wgspipeline/strelka:v0.0.1',
        'mutationseq': 'wgspipeline/mutationseq:v0.0.1',
        'vcftools': 'wgspipeline/vcftools:v0.0.1',
        'snpeff': 'wgspipeline/vcftools:v0.0.1',
        'titan': 'wgspipeline/titan:v0.0.1',
        'destruct': 'wgspipeline/destruct:v0.0.1',
        'lumpy': 'wgspipeline/lumpy:v0.0.1',
        'vizutils': 'wgspipeline/vizutils:v0.0.1',
        'fastqc': 'wgspipeline/fastqc:v0.0.1',
        'hmmcopy': 'wgspipeline/hmmcopy:v0.0.1',
    }

    singularity = {}

    return {'docker': docker_images, 'singularity': singularity}


def override_config(config, override):
    def update(d, u):
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    if not override:
        return config

    cfg = update(config, override)

    return cfg


def get_config(override):
    wgscontainers = containers()
    if override["cluster"] == "shahlab":
        config = shahlab_config(override["reference"], wgscontainers)
    elif override["cluster"] == "juno":
        config = juno_config(override["reference"], wgscontainers)
    elif override["cluster"] == "azure":
        config = azure_config(override["reference"], wgscontainers)
    else:
        raise Exception()

    override.pop('cluster')
    override.pop('reference')

    config = override_config(config, override)

    return config


def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.safe_dump(params, outputfile, default_flow_style=False)


def generate_pipeline_config(args):
    if args['which'] == 'generate_config':
        config_yaml = args['pipeline_config']
        config_yaml = os.path.abspath(config_yaml)
    else:
        config_yaml = "config.yaml"
        tmpdir = args.get("tmpdir", None)
        pipelinedir = args.get("pipelinedir", None)

        # use pypeliner tmpdir to store yaml
        if pipelinedir:
            config_yaml = os.path.join(pipelinedir, config_yaml)
        elif tmpdir:
            config_yaml = os.path.join(tmpdir, config_yaml)
        else:
            warnings.warn("no tmpdir specified, generating configs in working dir")
            config_yaml = os.path.join(os.getcwd(), config_yaml)

        config_yaml = helpers.get_incrementing_filename(config_yaml)

    params_override = {'cluster': 'azure', 'reference': 'grch37'}
    if args['config_override']:
        params_override.update(args["config_override"])

    helpers.makedirs(config_yaml, isfile=True)

    config = get_config(params_override)
    write_config(config, config_yaml)

    args["config_file"] = config_yaml

    return args
