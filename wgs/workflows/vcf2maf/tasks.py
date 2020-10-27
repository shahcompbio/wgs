import gzip
import os
import shutil

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
        tumour_id=None,
        normal_id=None,
        docker_image=None
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

    cmd = [
        'vcf2maf', vcf_unzipped, maf_output,
        os.path.join(reference, 'homo_sapiens', '99_GRCh37', 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'),
        os.path.join(reference, 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'),
        reference,
    ]

    if tumour_id:
        cmd.extend(['--tumor-id', tumour_id])
    if normal_id:
        cmd.extend(['--normal-id', normal_id])

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
