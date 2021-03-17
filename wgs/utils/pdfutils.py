import os

from PyPDF2 import PdfFileMerger
from fpdf import FPDF
from wgs.utils import helpers


def merge_titan_pngs(plots_dir, plot, num_clusters, chromosomes):
    imagelist = [os.path.join(plots_dir, "cluster_{}_chr{}.png".format(num_clusters, chrom))
                 for chrom in chromosomes]

    x, y, w, h = 0, 0, 200, 250

    pdf = FPDF()
    # imagelist is the list with all image filenames
    for image in imagelist:
        pdf.add_page()
        pdf.image(image, x, y, w, h)

    pdf.output(plot, "F")


def merge_pdfs(infiles, outfile):
    if isinstance(infiles, dict):
        infiles = infiles.values()

    merger = PdfFileMerger()

    for infile in infiles:
        # add it to list if not empty. skip empty files to avoid errors later
        if os.path.getsize(infile):
            merger.append(open(infile, 'rb'))

    helpers.makedirs(outfile, isfile=True)

    with open(outfile, 'wb') as fout:
        merger.write(fout)


def merge_jpgs(infiles, outfile):
    if isinstance(infiles, dict):
        infiles = infiles.values()

    x, y, w, h = 0, 0, 200, 250

    pdf = FPDF('P', 'mm', 'A4')
    # imagelist is the list with all image filenames
    for image in infiles:
        pdf.add_page()
        pdf.image(image, x, y, w, h)

    pdf.output(outfile, "F")
