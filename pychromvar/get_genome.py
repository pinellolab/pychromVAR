import os
import wget
import gzip
from tqdm import tqdm

def get_genome(genome:str="hg38", output_dir:str=None):
    """
    Download genome

    Args:
        genome (str, optional): 
            Which genome should be downloaded. Available options are: "hg19", "hg38", "mm9", "mm10".
            Defaults to "hg38".

        output_dir (str):
            Output directory. Default: current directory.
    """

    if not os.path.exists(output_dir):
        output_dir = os.getcwd()

    if genome == "hg19":
        _get_genome_hg19(output_dir=output_dir)


def _get_genome_hg19(output_dir):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
    chrom_list = ["chr" + str(e) for e in list(range(1, 23)) + ["X", "Y", "M"]]

    output_fname = os.path.join(output_dir, "genome_hg19.fa")
    with open(output_fname, "w") as f:
        for chrom in tqdm(chrom_list):
            gz_file_name = os.path.join(output_dir, chrom + ".fa.gz")
            wget.download(url + chrom + ".fa.gz", gz_file_name)

            gz_file = gzip.open(gz_file_name, "rb")
            f.write(gz_file.read().decode("utf-8"))
            gz_file.close()

            os.remove(gz_file_name)

