import pygenes


class InputArgsException(Exception):
    pass


class PygeneAnnotation(object):
    def __init__(self, infile, outfile, gtf_bin=None, gtf=None, is_contained=False):
        self.infile = infile
        self.outfile = outfile
        self.is_contained = is_contained

        if not gtf and not gtf_bin:
            raise InputArgsException('Requires either gtf or gtf_bin files')

        self.gene_models = pygenes.GeneModels()
        if gtf_bin:
            self.gene_models.load_binary(gtf_bin)
        else:
            self.gene_models.load_ensembl_gtf(gtf)

    def write_output(self):
        with open(self.infile, 'r') as titan_output, open(self.outfile, 'w') as writer:
            while True:
                line = titan_output.readline()
                if line[0] == '#':
                    continue
                header = line.rstrip()
                writer.write("%s\tPygenes(gene_id,gene_name;)\n" % header)
                break

            # do something to.. data lines
            for row in titan_output:
                row = row.rstrip()
                col = row.split('\t')

                chrom = col[1]
                start = int(col[2])
                end = int(col[3])

                if self.is_contained:
                    gene_ids = self.gene_models.find_contained_genes(chrom, start, end)
                else:
                    gene_ids = self.gene_models.find_overlapping_genes(chrom, start, end)

                pygenes_addition = ""
                for gene_id in gene_ids:
                    gene_name = self.gene_models.get_gene(gene_id).name

                    pygenes_addition += "%s,%s;" % (gene_id, gene_name)
                writer.write("%s\t%s\n" % (row, pygenes_addition))
