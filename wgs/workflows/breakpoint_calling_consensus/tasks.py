import pypeliner


def parse_destruct(infile, lumpy_data, output, low_map_filt_output, config, sample_id):
    '''
    Parse the input VCF file into a TSV file

    :param infile: temporary input VCF file
    :param output: path to the output TSV file
    '''

    cmd = ['vizutils_parse_destruct', '--infile', infile,
           '--pre_mappability_output', output,
           '--lumpy_data', lumpy_data,
           '--output', low_map_filt_output,
           '--case_id', sample_id,
           '--tumour_id', sample_id,
           '--normal_id', sample_id + 'N',
           ]

    for key, val in config.iteritems():
        if val is None:
            continue
        elif isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)
    pypeliner.commandline.execute(*cmd)


def parse_lumpy(infile, output, low_map_filt_output, config, sample_id):
    cmd = ['vizutils_parse_lumpy', '--infile', infile,
           '--pre_mappability_output', output,
           '--output', low_map_filt_output,
           '--case_id', sample_id]

    for key, val in config.iteritems():
        if val is None:
            continue
        elif isinstance(val, bool):
            if val:
                cmd.append('--{}'.format(key))
        else:
            cmd.append('--{}'.format(key))
            if isinstance(val, list):
                cmd.extend(val)
            else:
                cmd.append(val)
    pypeliner.commandline.execute(*cmd)
