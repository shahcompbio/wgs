import collections
import os
import warnings

import wgs
import yaml
from wgs.utils import helpers

from wgs.config import pipeline_config


def get_version():
    version = wgs.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


def containers():
    version = get_version()
    docker_images = {
        'bwa': 'bwa:v0.0.1',
        'samtools': 'samtools:v0.0.1',
        'picard': 'picard:v0.0.1',
        'wgs': 'wgs:v{}'.format(version),
        'strelka': 'strelka:v0.0.1',
        'mutationseq': 'mutationseq:v0.0.1',
        'vcftools': 'vcftools:v0.0.1',
        'snpeff': 'vcftools:v0.0.1',
        'titan': 'titan:v0.0.2',
        'destruct': 'destruct:v0.0.1',
        'lumpy': 'lumpy:v0.0.1',
        'vizutils': 'vizutils:v0.0.1',
        'fastqc': 'fastqc:v0.0.1',
        'hmmcopy': 'hmmcopy:v0.0.1',
    }

    return {'docker': docker_images}


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


def get_config(override, refdir):
    wgscontainers = containers()
    config = pipeline_config.pipeline_config(wgscontainers, refdir)
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

    helpers.makedirs(config_yaml, isfile=True)

    config = get_config(args['config_override'], args['refdir'])
    write_config(config, config_yaml)

    args["config_file"] = config_yaml

    return args
