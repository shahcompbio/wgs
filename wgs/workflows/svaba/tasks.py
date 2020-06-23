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
    shutil.copyfile(output_vcf + '.csi', outfile + '.csi')
    shutil.copyfile(output_vcf + '.tbi', outfile + '.tbi')


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
        docker_image=None
):
    cmd = svaba_cmd(tumor, normal, reference, tempdir, region=region, ncores=ncores, sample_id='sampleid')

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    get_vcfs(tempdir, 'sampleid.svaba.germline.indel.vcf', germline_indel)
    get_vcfs(tempdir, 'sampleid.svaba.germline.sv.vcf', germline_sv)
    get_vcfs(tempdir, 'sampleid.svaba.somatic.indel.vcf', somatic_indel)
    get_vcfs(tempdir, 'sampleid.svaba.somatic.sv.vcf', somatic_sv)

    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.germline.indel.vcf', unfiltered_germline_indel)
    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.germline.sv.vcf', unfiltered_germline_sv)
    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.somatic.indel.vcf', unfiltered_somatic_indel)
    get_vcfs(tempdir, 'sampleid.svaba.unfiltered.somatic.sv.vcf', unfiltered_somatic_sv)
