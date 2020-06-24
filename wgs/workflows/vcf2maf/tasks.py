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
    if vcf_file.endswith('.gz'):
        helpers.makedirs(tempdir)
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        gunzip_file(vcf_file, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file

    cmd = [
        'vcf2maf.pl', '--input-vcf', vcf_unzipped, '--output-maf', maf_output,
        '--vep-path', '/usr/local/bin',
        '--ref-fasta',
        os.path.join(reference, 'homo_sapiens', '99_GRCh37', 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'),
        '--filter-vcf', os.path.join(reference, 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'),
        '--vep-data', reference,
    ]

    if tumour_id:
        cmd.extend(['--tumor-id', tumour_id])
    if normal_id:
        cmd.extend(['--normal-id', normal_id])

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)
