import os

from fpdf import FPDF


def merge_pngs(plots_dir, plot, num_clusters, chromosomes):
    imagelist = [os.path.join(plots_dir, "cluster_{}_chr{}.png".format(num_clusters, chrom))
                 for chrom in chromosomes]

    x, y, w, h = 0, 0, 200, 250

    pdf = FPDF()
    # imagelist is the list with all image filenames
    for image in imagelist:
        pdf.add_page()
        pdf.image(image, x, y, w, h)

    pdf.output(plot, "F")
