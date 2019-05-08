class TransformVcfCounts(object):
    def __init__(self, infile, outfile, positions):
        self.infile = infile
        self.outfile = outfile
        self.positions = positions

    def read_ref_positions(self):
        if not self.args.positions_file:
            return

        ref_pos = set()
        freader = open(self.positions)
        for line in freader:
            line = line.strip().split(':')

            ref_pos.add(tuple(line))

        return ref_pos

    def main(self, ref_pos):

        with open(self.infile) as inputdata:
            with open(self.outfile, 'w') as outputdata:

                outputdata.write('chr\tposition\tref\trefCount\tNref\tNrefCount\n')

                for line in inputdata:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#"):
                        line = line.strip().split()
                        assert line[-1] == 'NORMAL', 'invalid vcf format'
                        assert line[-2] == 'TUMOUR', 'invalid vcf format'
                        continue
                    line = line.strip().split()
                    chrom = line[0]
                    pos = line[1]
                    ref = line[3]
                    assert line[8] == "RC:AC:NI:ND:DP:GT:PL", 'invalid vcf format'
                    tum_info = line[9].split(':')
                    tr = tum_info[0]
                    ta = tum_info[1]

                    if ref_pos and (chrom, pos) not in ref_pos:
                        continue

                    outstr = '\t'.join([chrom, pos, ref, tr, 'X', ta]) + '\n'
                    outputdata.write(outstr)
