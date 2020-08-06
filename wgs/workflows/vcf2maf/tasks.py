import gzip
import os
import shutil

import pypeliner
from wgs.utils import helpers

from wgs.utils import vcfutils

def gunzip_file(infile, outfile):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def run_vcf2maf(
        vcf_file,
        maf_output,
        tempdir,
        reference,
        vcf2maf_docker_image=None,
        vcftools_docker_image=None
):

    helpers.makedirs(tempdir)

    if vcf_file.endswith('.gz'):
        helpers.makedirs(tempdir)
        vcf_unzipped = os.path.join(tempdir, 'unzipped_vcf.vcf')
        gunzip_file(vcf_file, vcf_unzipped)
    else:
        vcf_unzipped = vcf_file

    sorted_file = os.path.join(tempdir, 'sorted.vcf')
    vcfutils.sort_vcf(vcf_unzipped, sorted_file, docker_image=vcftools_docker_image)

    cmd = [
        'vcf2maf.pl', '--input-vcf', sorted_file, '--output-maf', maf_output,
        '--vep-path', '/usr/local/bin',
        '--ref-fasta',
        os.path.join(reference, 'homo_sapiens', '99_GRCh37', 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'),
        '--filter-vcf', os.path.join(reference, 'ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz'),
        '--vep-data', reference,
    ]

    cmd = ' '.join(map(str,cmd))
    pypeliner.commandline.execute(*cmd, docker_image=vcf2maf_docker_image)

def update_maf_counts(input_maf, counts_file, output_maf, sample_id):
    counts = {}
    with open(counts_file) as infile:
        for line in infile:
            line = line.strip().split()
            chrom, pos, id, ta, tr, td, na, nr, nd = line
            counts[(chrom, pos, id)] = (ta, tr, td, na, nr, nd)

    with open(input_maf) as infile, open(output_maf, 'wt') as outfile:
        header = infile.readline()
        assert header.startswith('#')
        outfile.write(header)

        header = infile.readline()
        outfile.write(header)

        header = {v: i for i, v in enumerate(header.strip().split('\t'))}

        t_dp = header['t_depth']
        t_ref = header['t_alt_count']
        t_alt = header['t_ref_count']
        n_dp = header['n_depth']
        n_ref = header['n_alt_count']
        n_alt = header['n_ref_count']

        tum_sample_barcode = header['Tumor_Sample_Barcode']
        norm_sample_barcode = header['Matched_Norm_Sample_Barcode']

        chrom = header['Chromosome']
        pos = header['vcf_pos']
        vcfid = header['vcf_id']

        for line in infile:
            line_split = line.strip().split('\t')

            if (line_split[chrom], line_split[pos], line_split[vcfid]) in counts:
                ta, tr, td, na, nr, nd = counts[(line_split[chrom], line_split[pos], line_split[vcfid])]

                line_split[t_dp] = td
                line_split[t_ref] = tr
                line_split[t_alt] = ta

                line_split[n_dp] = nd
                line_split[n_ref] = nr
                line_split[n_alt] = na

                line_split[tum_sample_barcode] = sample_id
                line_split[norm_sample_barcode] = sample_id+'_N'

                line = '\t'.join(line_split) + '\n'

            outfile.write(line)
