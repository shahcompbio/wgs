'''
Created on Feb 19, 2018

@author: dgrewal
'''
import collections
import errno
import gzip
import logging
import multiprocessing
import os
import re
import shutil
import tarfile
from multiprocessing.pool import ThreadPool
from subprocess import Popen, PIPE

import pypeliner
import yaml

import wgs


class GetFileHandle(object):
    def __init__(self, filename, mode='rt'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        if self.get_file_format(self.filename) in ["csv", 'plain-text']:
            self.handle = open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "gzip":
            self.handle = gzip.open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "h5":
            self.handle = pd.HDFStore(self.filename, self.mode)
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    @property
    def handler(self):
        return self.__enter__()

    def get_file_format(self, filepath):
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        else:
            logging.getLogger("wgs.helpers").warning(
                "Couldn't detect output format. extension {}".format(ext)
            )
            return "plain-text"


def get_default_ctx(memory=4, walltime='04:00', ncpus=1, disk=8):
    ctx = {
        'mem': memory,
        'walltime': walltime,
        'mem_num_retry': 3,
        'mem_retry_factor': 2,
        'walltime_num_retry': 5,
        'walltime_retry_increment': '24:00',
        'ncpus': ncpus,
        'disk': disk,
    }

    return ctx


def get_fastqs(inputs):
    fq1 = {}
    fq2 = {}

    for lane in inputs['fastqs']:
        print(lane)
        print(inputs[lane])
        fq1[lane] = inputs[lane]['fastq1']
        fq2[lane] = inputs[lane]['fastq2']

    return fq1, fq2


def build_shell_script(command, tag, tempdir):
    outfile = os.path.join(tempdir, "{}.sh".format(tag))
    with open(outfile, 'w') as scriptfile:
        scriptfile.write("#!/bin/bash\n")
        if isinstance(command, list) or isinstance(command, tuple):
            command = ' '.join(map(str, command)) + '\n'
        scriptfile.write(command)
    return outfile


def run_in_gnu_parallel(commands, tempdir, ncores=None):
    makedirs(tempdir)

    scriptfiles = []

    for tag, command in enumerate(commands):
        scriptfiles.append(build_shell_script(command, tag, tempdir))

    parallel_outfile = os.path.join(tempdir, "commands.txt")
    with open(parallel_outfile, 'w') as outfile:
        for scriptfile in scriptfiles:
            outfile.write("sh {}\n".format(scriptfile))

    if not ncores:
        ncores = multiprocessing.cpu_count()

    gnu_parallel_cmd = ['parallel', '--jobs', ncores, '<', parallel_outfile]
    pypeliner.commandline.execute(*gnu_parallel_cmd)


def get_values_from_input(yamldata, key):
    values = {str(sample): yamldata[sample].get(key) for sample in yamldata}
    return values


def get_coltype_reference():
    coltypes = {
        'estimated_library_size': 'int', 'total_mapped_reads': 'int',
        'total_reads_hmmcopy': 'float', 'cor_map': 'float',
        'cv_neutral_state': 'float', 'MBRSM_dispersion': 'float',
        'map': 'float', 'mad_hmmcopy': 'float', 'copy': 'float',
        'modal_curve': 'float', 'sample_well': 'str', 'true_multiplier': 'float',
        'reads': 'int', 'jira_id': 'str', 'gc': 'float', 'integer_copy_number': 'int',
        'breakpoints': 'int', 'total_duplicate_reads': 'int',
        'quality': 'float', 'cv_hmmcopy': 'float', 'empty_bins_hmmcopy': 'float',
        'paired_mapped_reads': 'int', 'total_reads': 'int', 'end': 'int', 'width': 'int',
        'MBRSI_dispersion_non_integerness': 'float', 'total_properly_paired': 'int',
        'chr': 'str', 'sample_type': 'str', 'mean_insert_size': 'float', 'start': 'int',
        'state': 'int', 'valid': 'bool', 'coverage_breadth': 'float', 'empty_bins_hmmcopy_chrY': 'int',
        'too_even': 'bool', 'unpaired_duplicate_reads': 'int', 'unpaired_mapped_reads': 'int',
        'unmapped_reads': 'int', 'mad_chr19': 'float', 'cell_id': 'str', 'cell_call': 'str',
        'coverage_depth': 'float', 'median_insert_size': 'float', 'modal_quantile': 'float',
        'sample_plate': 'str', 'mean_state_mads': 'float', 'ideal': 'bool',
        'experimental_condition': 'str', 'mean_copy': 'float', 'mean_hmmcopy_reads_per_bin': 'float',
        'multiplier': 'int', 'percent_duplicate_reads': 'float', 'i7_barcode': 'str',
        'total_halfiness': 'float', 'std_hmmcopy_reads_per_bin': 'float',
        'standard_deviation_insert_size': 'float', 'mean_state_vars': 'float',
        'all_heatmap_order': 'int', 'scaled_halfiness': 'float', 'cor_gc': 'float',
        'median': 'float', 'state_mode': 'int', 'paired_duplicate_reads': 'int',
        'median_hmmcopy_reads_per_bin': 'float', 'mad_neutral_state': 'float',
        'autocorrelation_hmmcopy': 'float', 'mad_autosomes': 'float', 'i5_barcode': 'str',
        'loglikehood': 'float', 'MSRSI_non_integerness': 'float'}

    ignore_cols = set(range(9))

    return coltypes, ignore_cols


def resolve_template(regions, template, format_key):
    outputs = {v: template.format(**{format_key: v}) for v in regions}
    return outputs


def get_fastq_files(input_yaml):
    data = load_yaml(input_yaml)

    items = {}
    for lane, laneinfo in data["fastqs"].items():
        items[lane] = {}
        items[lane]['fastq_1'] = format_file_yaml(laneinfo['fastq_1'])
        items[lane]['fastq_2'] = format_file_yaml(laneinfo['fastq_2'])

    return items


def format_file_yaml(filepath):
    ext = os.path.splitext(filepath)

    if ext[1] == ".gz":
        ext = os.path.splitext(ext[0])

    mapping = {'.bam': 'bam', '.pdf': 'PDF',
               '.fastq': 'fastq', '.h5': 'H5',
               '.tar': 'tar'}

    return {'filename': filepath, 'type': mapping[ext[1]]}


def write_to_yaml(outfile, data):
    with open(outfile, 'w') as output:
        yaml.safe_dump(data, output, default_flow_style=False)


def eval_expr(val, operation, threshold):
    if operation == "gt":
        if val > threshold:
            return True
    elif operation == 'ge':
        if val >= threshold:
            return True
    elif operation == 'lt':
        if val < threshold:
            return True
    elif operation == 'le':
        if val <= threshold:
            return True
    elif operation == 'eq':
        if val == threshold:
            return True
    elif operation == 'ne':
        if not val == threshold:
            return True
    else:
        raise Exception("unknown operator type: {}".format(operation))

    return False


def get_incrementing_filename(path):
    """
    avoid overwriting files. if path exists then return path
    otherwise generate a path that doesnt exist.
    """

    if not os.path.exists(path):
        return path

    i = 0
    while os.path.exists("{}.{}".format(path, i)):
        i += 1

    return "{}.{}".format(path, i)


def run_in_parallel(worker, args, ncores=None):
    def args_unpack(worker, args):
        return worker(*args)

    count = multiprocessing.cpu_count()

    if ncores:
        count = min(ncores, count)

    pool = ThreadPool(processes=count)

    tasks = []

    for arg in args:
        task = pool.apply_async(args_unpack,
                                args=(worker, arg),
                                )
        tasks.append(task)

    pool.close()
    pool.join()

    [task.get() for task in tasks]

    pool.terminate()
    del pool


def run_cmd(cmd, output=None):
    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(
                cmdout,
                cmderr))

    if output:
        stdout.close()


def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.safe_load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data


def symlink(actual_file, symlink):
    if not os.path.exists(symlink):
        os.symlink(actual_file, symlink)


def copy_file(infile, output):
    shutil.copy(infile, output)


def get_instrument_info(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.items():

            if "sequencing_instrument" not in paths:
                raise Exception(
                    "instrument key missing in cell: {}".format(cell))
            seqinfo[(cell, lane)] = paths["sequencing_instrument"]

    return seqinfo


def get_center_info(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "fastqs" in data[
            cell], "couldnt extract fastq file paths from yaml input for cell: {}".format(cell)

    seqinfo = dict()
    for cell in data.keys():
        fastqs = data[cell]["fastqs"]

        for lane, paths in fastqs.items():

            if "sequencing_center" not in paths:
                raise Exception(
                    "sequencing_center key missing in cell: {}".format(cell))
            seqinfo[(cell, lane)] = paths["sequencing_center"]

    return seqinfo


def get_samples(fastqs_file):
    data = load_yaml(fastqs_file)

    return list(data.keys())


def get_bams(fastqs_file):
    data = load_yaml(fastqs_file)

    for cell in data.keys():
        assert "bam" in data[
            cell], "couldnt extract bam file paths from yaml input for cell: {}".format(cell)

    bam_filenames = {cell: data[cell]["bam"] for cell in data.keys()}
    bai_filenames = {cell: data[cell]["bam"] + ".bai" for cell in data.keys()}

    return bam_filenames, bai_filenames


def load_config(args):
    return load_yaml(args["config_file"])


def makedirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def rmdirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)

    if not os.path.exists(directory):
        return

    try:
        shutil.rmtree(directory)
    except OSError:
        if os.path.exists(directory):
            raise


def add_extensions(filepaths):
    paths_extensions = []
    for filepath in filepaths:
        paths_extensions.append(filepath)

        if filepath.endswith('.csv.gz') or filepath.endswith('csv'):
            paths_extensions.append(filepath + '.yaml')
        elif filepath.endswith('.vcf.gz'):
            paths_extensions.append(filepath + '.csi')
            paths_extensions.append(filepath + '.tbi')
        elif filepath.endswith('.bam'):
            paths_extensions.append(filepath + '.bai')

    return paths_extensions


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def make_tar_from_files(output_filename, input_files, tempdir):
    makedirs(tempdir)

    for infile in input_files:
        if isinstance(infile, str):
            shutil.copyfile(infile, os.path.join(tempdir, os.path.basename(infile)))
        elif isinstance(infile, collections.Mapping):
            for key, filepath in infile.items():
                if not isinstance(key, str):
                    key = ' '.join(key)
                shutil.copyfile(filepath, os.path.join(tempdir, '{}_{}'.format(key, os.path.basename(filepath))))

    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(tempdir, arcname=os.path.basename(tempdir))


def generate_meta_yaml_file(
        metadata_file,
        filepaths=None,
        metadata=None,
        root_dir=None
):
    if not root_dir:
        final_paths = filepaths
    else:
        final_paths = []
        for filepath in filepaths:
            if not filepath.startswith(root_dir):
                error = 'file {} does not have {} in path'.format(
                    filepath, root_dir
                )
                raise Exception(error)

            filepath = os.path.relpath(filepath, root_dir)
            final_paths.append(filepath)

    final_paths = add_extensions(final_paths)

    metadata = {
        'filenames': final_paths,
        'meta': metadata,
    }

    write_to_yaml(metadata_file, metadata)


def expand_list(list, expanders, to_replace):
    '''
    for each str in list,
    creates expander items,
    each marked with expander
    if {to_replace} in list,
    to_replace is replaced
    with expander items

    :param list: list of strs
    :param expanders: list of markers for list
    :param to_replace: item to replace in each of
    the items in list
    '''
    outlist = []
    replace_list = {to_replace: ''}
    for item in list:
        if isinstance(item, dict):
            outlist.extend(item.values())
        if to_replace in item:
            for expander in expanders:
                replace_list[to_replace] = expander
                final_name = item.format(**replace_list)
                outlist.append(final_name)
    return outlist


def get_version():
    '''
    gets the version of wgs
    '''
    version = wgs.__version__
    # strip setuptools metadata
    version = version.split("+")[0]
    return version


class InputException(Exception):
    pass


def generate_and_upload_metadata(
        command, root_dir, filepaths, output, template=None,
        input_yaml_data=None, input_yaml=None, metadata={}, type=None
):
    if not metadata:
        metadata = {}

    if isinstance(filepaths, dict):
        filepaths = filepaths.values()
    filepaths = list(filepaths)

    metadata['command'] = ' '.join(command)
    metadata['version'] = get_version()

    if type:
        metadata['type'] = type

    if template:
        assert len(template) == 3
        instances, template_path, instance_key = template
        assert re.match('.*\{.*\}.*', template_path)
        template_path = os.path.relpath(template_path, root_dir)
        metadata['bams'] = {}
        metadata['bams']['template'] = template_path
        instances = [{instance_key: instance} for instance in instances]
        metadata['bams']['instances'] = instances

    if input_yaml_data:
        if not input_yaml:
            raise InputException("missing yaml file to write to")
        with open(input_yaml, 'wt') as yaml_writer:
            yaml.safe_dump(input_yaml_data, yaml_writer)

        if not input_yaml.startswith(root_dir) and root_dir in input_yaml:
            input_yaml = input_yaml[input_yaml.index(root_dir):]
        if input_yaml.endswith('.tmp'):
            input_yaml = input_yaml[:-4]

        metadata['input_yaml'] = os.path.relpath(input_yaml, root_dir)
        filepaths.append(input_yaml)

    generate_meta_yaml_file(
        output, filepaths=filepaths, metadata=metadata, root_dir=root_dir
    )


def load_qc_input_yaml_flat(path):
    data = {}
    yaml = load_yaml(path)

    for cohort, cohort_data in yaml.items():
        for sample, sample_data in cohort_data.items():
            data[(cohort, sample)] = {data_label: data
                                      for data_label, data in sample_data.items()
                                      }
    return data
