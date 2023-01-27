pychromVAR
==============================================================

pychromVAR is a python package for inferring transcription factor binding variability from 
scATAC-seq data by implmenting the algorithm proposed in [chromVAR](https://github.com/GreenleafLab/chromVAR). 
It is built on [anndata](https://anndata.readthedocs.io/en/latest/) and 
[mudata](https://mudata.readthedocs.io/en/latest/) therefore can work seamlessly with [Scanpy](https://scanpy.readthedocs.io/en/stable/) and [Muon](https://muon.readthedocs.io/en/latest/) pipeline.

For more methdological detials, please refer to the original [paper](https://www.nature.com/articles/nmeth.4401).

.. toctree::
   :caption: MAIN
   :maxdepth: 2
   :hidden:

   installation

.. toctree::
   :caption: NOTEBOOKS
   :maxdepth: 2
   :hidden:

   notebooks/run_chromVAR