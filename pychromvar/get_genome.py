import os
import wget
import gzip
from tqdm import tqdm

base_url = "https://ftp.ebi.ac.uk/pub/databases/gencode"

genome_url = {
    "hg19": f"{base_url}/Gencode_human/release_42/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz",
    "hg38": f"{base_url}/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz",
    "mm10": f"{base_url}/Gencode_mouse/release_M1/NCBIM37.genome.fa.gz",
    "mm39": f"{base_url}/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz"
}

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
    assert genome in ["hg19", "hg38", "mm10", "mm39"], f"Cannot find {genome}!"

    if not os.path.exists(output_dir):
        output_dir = os.getcwd()

    output_fname = os.path.join(output_dir, f"{genome}.fa")
    if os.path.exists(output_fname):
        os.remove(output_fname)

    gz_fname = os.path.join(output_dir, f"{genome}.fa.gz")
    wget.download(genome_url[genome], gz_fname, bar=bar_thermometer)

    output_fname = os.path.join(output_dir, f"{genome}.fa")
    with open(output_fname, "w") as f:
        gz_file = gzip.open(gz_fname, "rb")
        f.write(gz_file.read().decode("utf-8"))
        gz_file.close()

    os.remove(gz_fname)
