import pypeliner


def parse_destruct(infile, lumpy_data, output, low_map_filt_output, config):
    '''
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    '''

    cmd = ['vizutils_parse_destruct', '--infile', infile,
           '--pre_mappability_output', output,
           '--lumpy_data', lumpy_data,
           '--output', low_map_filt_output]

    for key, val in config.iteritems():
        if val is None:
            continue
        elif isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            if isinstance(val, list):
                val = ' '.join(val)

            cmd.extend(['--{}'.format(key), val])

    pypeliner.commandline.execute(*cmd)


def parse_lumpy(infile, output, low_map_filt_output, config):
    cmd = ['vizutils_parse_lumpy', '--infile', infile,
           '--pre_mappability_output', output,
           '--output', low_map_filt_output,]

    for key, val in config.iteritems():
        if val is None:
            continue
        elif isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            if isinstance(val, list):
                val = ' '.join(val)
            cmd.extend(['--{}'.format(key), val])

    pypeliner.commandline.execute(*cmd)

