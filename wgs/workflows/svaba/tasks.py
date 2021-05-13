import os
import shutil

import pypeliner
from wgs.utils import helpers


def svaba_cmd(tumor, normal, reference, tempdir, region=None, ncores=None, sample_id='sample'):
    helpers.makedirs(tempdir)
    tempdir = os.path.join(tempdir, sample_id)

    cmd = [
        'svaba', 'run', '-t', tumor, '-n', normal,
        '-G', reference, '-z', '-a', tempdir
    ]

    if region:
        cmd += ['-k', region]

    if ncores:
        cmd += ['-p', ncores]

    return cmd


def get_vcfs(tempdir, filename, outfile):
    output_vcf = os.path.join(tempdir, filename)
    shutil.copyfile(output_vcf, outfile)


def run_svaba(
        tumor,
        normal,
        germline_indel,
        germline_sv,
        somatic_indel,
        somatic_sv,
        unfiltered_germline_indel,
        unfiltered_germline_sv,
        unfiltered_somatic_indel,
        unfiltered_somatic_sv,
        reference,
        tempdir,
        region=None,
        ncores=None,
):
    cmd = svaba_cmd(tumor, normal, reference, tempdir, region=region, ncores=ncores, sample_id='sampleid')

    pypeliner.commandline.execute(*cmd)

    get_vcfs(tempdir, 'sampleid.svaba.germline.indel.vcf.gz', germline_indel)
    get_vcfs(tempdir, 'sampleid.svaba.germline.sv.vcf.gz', germline_sv)
    get_vcfs(tempdir, 'sampleid.svaba.somatic.indel.vcf.gz', somatic_indel)
    get_vcfs(tempdir, 'sampleid.svaba.somatic.sv.vcf.gz', somatic_sv)

    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.germline.indel.vcf.gz', unfiltered_germline_indel)
    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.germline.sv.vcf.gz', unfiltered_germline_sv)
    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.somatic.indel.vcf.gz', unfiltered_somatic_indel)
    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.somatic.sv.vcf.gz', unfiltered_somatic_sv)
