[![Stars](https://img.shields.io/github/stars/lzj1769/pychromVAR?logo=GitHub&color=yellow)](https://github.com/lzj1769/pychromVAR/stargazers)
[![PyPI](https://img.shields.io/pypi/v/pychromvar?logo=PyPI)](https://pypi.org/project/pychromvar/)
[![Docs](https://readthedocs.org/projects/pychromvar/badge/?version=latest)](https://pychromvar.readthedocs.io)

# pychromVAR 

pychromVAR is a python package for inferring transcription factor binding variability from scATAC-seq data by implmenting the algorithm proposed in [chromVAR](https://github.com/GreenleafLab/chromVAR). It is built on [anndata](https://anndata.readthedocs.io/en/latest/) and [mudata](https://mudata.readthedocs.io/en/latest/) therefore can work seamlessly with [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) and [Muon](https://gtca.github.io/muon/) pipeline. 

For more methodological details, please refer to the [paper](https://www.nature.com/articles/nmeth.4401). 

# Installation

The quickest and easiest way to get pychromvar is to to use pip:

```shell
pip install pychromvar
```

# Tutorial

You can find [here](https://pychromvar.readthedocs.io/en/latest/notebooks/multimodal_pbmc_3k.html) a tutorial how to use pychromVAR combined with Muon to analysis multimodal single-cell PBMC data.


