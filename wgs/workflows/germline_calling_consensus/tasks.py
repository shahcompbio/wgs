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
