import gzip
import os
import shutil

import pandas as pd
import pypeliner
from wgs.utils import helpers
import vcf
import itertools

def gunzip_file(infile, outfile):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def run_vcf2maf(
        vcf_file,
        maf_output,
        tempdir,
        reference,
        vep_fasta_suffix,
        vep_ncbi_build,
        vep_cache_version,
):
    if os.path.exists(tempdir):
        helpers.rmdirs(tempdir)

    helpers.makedirs(tempdir)

    input_vcf = os.path.join(tempdir, os.path.basename(vcf_file))
    shutil.copyfile(vcf_file, input_vcf)

    if vcf_file.endswith('.gz'):
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        gunzip_file(input_vcf, vcf_unzipped)
    else:
        vcf_unzipped = input_vcf

    assert vcf_unzipped.endswith('.vcf')
    vcf_unzipped_vep = vcf_unzipped[:-4]
    vcf_unzipped_vep = vcf_unzipped_vep+'.vep.vcf'

    if os.path.exists(vcf_unzipped_vep):
        os.remove(vcf_unzipped_vep)


    cmd = [
        'vcf2maf', vcf_unzipped, maf_output,
        os.path.join(reference, vep_fasta_suffix),
        reference, vep_ncbi_build, vep_cache_version
    ]

    pypeliner.commandline.execute(*cmd)


def update_ids(infile, tumour_id, normal_id, output):
    with open(infile) as infile_read:
        maf_header = infile_read.readline()
    assert maf_header.startswith('#version 2.4')

    df = pd.read_csv(infile, dtype='str', skiprows=1, sep='\t')

    if len(df)>1:
        assert len(df['Tumor_Sample_Barcode'].unique()) == 1
        assert df['Tumor_Sample_Barcode'].unique()[0] == 'TUMOR'

        assert len(df['Matched_Norm_Sample_Barcode'].unique()) == 1
        assert df['Matched_Norm_Sample_Barcode'].unique()[0] == 'NORMAL'

        df['Matched_Norm_Sample_Barcode'] = normal_id

        # for germlines tumour will be none
        if tumour_id is None:
            tumour_id = 'NA'
        df['Tumor_Sample_Barcode'] = tumour_id

    with open(output, 'wt') as outfile:
        outfile.write(maf_header)

        df.to_csv(outfile, sep='\t', index=False)


def split_vcf(in_file, out_file_callback, lines_per_file):
    """ Split a VCF file into smaller files.

    :param in_file: Path of VCF file to split.

    :param out_file_callback: Callback function which supplies file name given index of split.

    :param lines_per_file: Maximum number of lines to be written per file.

     """

    def line_group(_, line_idx=itertools.count()):
        return int(next(line_idx) / lines_per_file)

    reader = vcf.Reader(filename=in_file)

    for file_idx, records in itertools.groupby(reader, key=line_group):
        file_name = out_file_callback(file_idx)

        with open(file_name, 'wt') as out_fh:
            writer = vcf.Writer(out_fh, reader)

            for record in records:
                writer.write_record(record)

            writer.close()


def merge_mafs(maf_files, output):

    if isinstance(maf_files, dict):
        maf_files = list(maf_files.values())

    with helpers.GetFileHandle(output, 'wt') as maf_writer:

        with helpers.GetFileHandle(maf_files[0]) as header_read:
            header = header_read.readline()
            assert header.startswith('#version 2.4')
            maf_writer.write(header)

            header = header_read.readline()
            assert header.startswith('Hugo_Symbol')
            maf_writer.write(header)

        for filepath in maf_files:
            with helpers.GetFileHandle(filepath, 'rt') as maf_reader:
                for line in maf_reader:
                    if line.startswith('Hugo_Symbol') or line.startswith('#'):
                        continue
                    maf_writer.write(line)
