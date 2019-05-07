import pypeliner


def parse_museq(infile, output, low_map_filt_output, config, docker_image=None):
    '''
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    '''

    cmd = ['vizutils_parse_museq', '--infile', infile,
           '--pre_mappability_output', output,
           '--output', low_map_filt_output]

    for key, val in config.iteritems():
        if isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def parse_strelka(infile, output, low_map_filt_output, config, docker_image=None):
    cmd = ['vizutils_parse_strelka', '--infile', infile,
           '--pre_mappability_output', output,
           '--output', low_map_filt_output, ]

    for key, val in config.iteritems():
        if isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)
    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def merge_overlap(infiles, outfile, docker_image=None):
    cmd = ['vizutils_merge', '--output', outfile, 'merge_overlap',
           '--input1', infiles[0], '--input2', infiles[1],
           '--suffix1', 'strelka', '--key_cols', 'case_id',
           'chromosome', 'start', 'stop', 'ref', 'alt']

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)


def concatenate(infiles, outfile, docker_image=None):
    cmd = ['vizutils_merge', '--output', outfile, 'concatenate',
           '--input', ]

    cmd += infiles

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
