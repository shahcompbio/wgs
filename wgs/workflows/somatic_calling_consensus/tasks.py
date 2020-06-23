def update_maf_counts(input_maf, counts_file, output_maf):
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

                line = '\t'.join(line_split) + '\n'

            outfile.write(line)
