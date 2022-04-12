import gzip
import os
import shutil

import pandas as pd
import pypeliner
from wgs.utils import helpers


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
