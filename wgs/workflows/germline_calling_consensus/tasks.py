from wgs.utils import helpers


def update_maf_counts(input_maf, counts_file, output_maf):
    counts = {}
    with open(counts_file) as infile:
        for line in infile:
            line = line.strip().split()
            chrom, pos, id, na, nr, nd = line
            counts[(chrom, pos, id)] = (na, nr, nd)

    with open(input_maf) as infile, open(output_maf, 'wt') as outfile:
        header = infile.readline()
        assert header.startswith('#')
        outfile.write(header)

        header = infile.readline()
        outfile.write(header)

        header = {v: i for i, v in enumerate(header.strip().split('\t'))}
        n_dp = header['n_depth']
        n_ref = header['n_alt_count']
        n_alt = header['n_ref_count']

        chrom = header['Chromosome']
        pos = header['vcf_pos']
        vcfid = header['vcf_id']

        for line in infile:
            line = line.strip().split('\t')

            try:
                na, nr, nd = counts[(line[chrom], line[pos], line[vcfid])]
            except:
                print(line)
                print(chrom, pos)
                print(line[chrom])
                print(line[pos])
                raise

            line[n_dp] = nd
            line[n_ref] = nr
            line[n_alt] = na

            line = '\t'.join(line) + '\n'

            outfile.write(line)


def split_vcf_by_chr(vcf_file, chromosome, output):
    with helpers.GetFileHandle(vcf_file, 'rt') as vcf_reader, \
            helpers.GetFileHandle(output, 'wt') as vcf_writer:
        for line in vcf_reader:
            if line.startswith('#') or line.startswith(chromosome):
                vcf_writer.write(line)


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
