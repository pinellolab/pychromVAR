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
    assert genome in ["hg19", "hg38", "mm9", "mm10"], f"Cannot find {genome}!"

    if not os.path.exists(output_dir):
        output_dir = os.getcwd()

    output_fname = os.path.join(output_dir, f"{genome}.fa")
    if os.path.exists(output_fname):
        os.remove(output_fname)

    if genome == "hg19":
        _get_genome_hg19(output_dir=output_dir)
    elif genome == "hg38":
        _get_genome_hg38(output_dir=output_dir)
    elif genome == "mm9":
        _get_genome_mm9(output_dir=output_dir)
    elif genome == "mm10":
        _get_genome_mm10(output_dir=output_dir)


def _get_genome_hg19(output_dir):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
    chrom_list = ["chr" + str(e) for e in list(range(1, 23)) + ["X", "Y", "M"]]

    output_fname = os.path.join(output_dir, "hg19.fa")
    with open(output_fname, "w") as f:
        for chrom in tqdm(chrom_list):
            gz_file_name = os.path.join(output_dir, chrom + ".fa.gz")
            wget.download(url + chrom + ".fa.gz", gz_file_name, bar=None)
            gz_file = gzip.open(gz_file_name, "rb")
            f.write(gz_file.read().decode("utf-8"))

            gz_file.close()
            os.remove(gz_file_name)

def _get_genome_hg38(output_dir):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/"
    chrom_list = ["chr" + str(e) for e in list(range(1, 23)) + ["X", "Y", "M"]]

    output_fname = os.path.join(output_dir, "hg38.fa")
    with open(output_fname, "w") as f:
        for chrom in tqdm(chrom_list):
            gz_file_name = os.path.join(output_dir, chrom + ".fa.gz")
            wget.download(url + chrom + ".fa.gz", gz_file_name, bar=None)
            gz_file = gzip.open(gz_file_name, "rb")
            f.write(gz_file.read().decode("utf-8"))

            gz_file.close()
            os.remove(gz_file_name)

def _get_genome_mm9(output_dir):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/"
    chrom_list = ["chr" + str(e) for e in list(range(1, 29)) + ["X", "Y", "M"]]

    output_fname = os.path.join(output_dir, "mm9.fa")
    with open(output_fname, "w") as f:
        for chrom in tqdm(chrom_list):
            gz_file_name = os.path.join(output_dir, chrom + ".fa.gz")
            wget.download(url + chrom + ".fa.gz", gz_file_name, bar=None)
            gz_file = gzip.open(gz_file_name, "rb")
            f.write(gz_file.read().decode("utf-8"))

            gz_file.close()
            os.remove(gz_file_name)

def _get_genome_mm10(output_dir):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/"
    chrom_list = ["chr" + str(e) for e in list(range(1, 29)) + ["X", "Y", "M"]]

    output_fname = os.path.join(output_dir, "mm10.fa")
    with open(output_fname, "w") as f:
        for chrom in tqdm(chrom_list):
            gz_file_name = os.path.join(output_dir, chrom + ".fa.gz")
            wget.download(url + chrom + ".fa.gz", gz_file_name, bar=None)
            gz_file = gzip.open(gz_file_name, "rb")
            f.write(gz_file.read().decode("utf-8"))

            gz_file.close()
            os.remove(gz_file_name)