# BenchmarkST

We present a comprehensive evaluation of the performance of state-of-the-art methods for cell/spot clustering, alignment, and integration on spatial transcriptomic datasets. Our goal is to provide the research field with a fair evaluation and comparison of all available methods to facilitate the users’ choice as well as shed more light on further developments to focus on the challenging aspects of ST based tasks.

## Table of Content 
- [Online Documentation](#online-documentationtutorials)
- [Installation](#installation)
- [Repository description](#Repository-description)
- [General Usage](#General-Usage)
  -  [ADEPT](#ADEPT)
  -  [GraphST](#GraphST)
  -  [conST](#conST)
  -  [conGI](#conGI)
  -  [SpatialPCA](#SpatialPCA)
  -  [DR-SC](#DR-SC)
  -  [STAGATE](#STAGATE)
  -  [CCST](#CCST)
  -  [SEDR](#SEDR)
  -  [SpaGCN](#SpaGCN)
- [Citation](#Citation)

## Online documentation&tutorials

We hope the online readthedocs.io version of doc of our project BenchmarkST would make things easier and more organized. Feel free to navigate to it through the following link.

Online Doc: [Here](https://benchmarkst-reproducibility.readthedocs.io/en/latest/index.html)


## Installation

We strongly suggest to follow the installation instructions given by each method. We will also provide the .yaml file of our debugging environment in each subfolder.

## Repository description
We provide all the jupyter notebooks in each subfolder to run and evaluate all algorithms for all datasets with ground truth annotation, and to reproduce the results introduced in our paper.


<!-- 1. '${algorithm_name}' folders contain mostly the original repo's content and code, but we also adapt a bit so that it could be run on the new datasets that have not been tested in their published paper.
2. 'analysis' folder constains script for reproducing the figures in our paper.
3. -------
4. ------- -->

## General Usage

Please follow the readme or jupyter notebook instructions in each subfolder for further details.

### original Github repository links and our script for Clustering task:

#### BANKSY

Github link:https://github.com/prabhakarlab/Banksy

Our tutorial: [run_BANKSY]()

#### ADEPT

Github link: https://github.com/maiziezhoulab/ADEPT

Our tutorial: [run_ADEPT](https://github.com/maiziezhoulab/BenchmarkST/blob/main/ADEPT/run_ADEPT/run_ADEPT.ipynb)

#### GraphST

Github link: https://github.com/JinmiaoChenLab/GraphST

<!-- installation: https://deepst-tutorials.readthedocs.io/en/latest/Installation.html -->

Our tutorial: [run_GraphST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/GraphST/run_GraphST/tutorial.ipynb)

#### SpaceFLow

Github link: https://github.com/hongleir/SpaceFlow

Our tutorial: [run_SpaceFLow]()


#### conST

Github link: https://github.com/ys-zong/conST

Our tutorial: [run_conST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/conST/run_conST.ipynb)

#### ConGI

Github link: https://github.com/biomed-AI/ConGI

Our tutorial: [run_conGI](https://github.com/maiziezhoulab/BenchmarkST/blob/main/ConGI/run_ConGI/run_conGI.ipynb)

#### SpatialPCA

Github link: https://github.com/shangll123/SpatialPCA

Our tutorial: [run_SpatialPCA]()

#### DR-SC

Github link: https://github.com/feiyoung/DR.SC

Our tutorial: [run_DRSC]()

#### STAGATE

Github link: https://github.com/QIFEIDKN/STAGATE_pyG

<!-- installation: https://stagate.readthedocs.io/en/latest/Installation_pyG.html -->

Our tutorial: [run_STAGATE](https://github.com/maiziezhoulab/BenchmarkST/blob/main/STAGATE_pyG/run_STAGATE/run_STAGATE.ipynb)

#### CCST

Github link: https://github.com/xiaoyeye/CCST

Our tutorial: [run_CCST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/CCST/run_CCST.ipynb)

#### SEDR

Github link: https://github.com/JinmiaoChenLab/SEDR

Our tutorial: [run_SEDR](https://github.com/maiziezhoulab/BenchmarkST/blob/main/SEDR/run_SEDR.ipynb)

#### SpaGCN

Github link: https://github.com/jianhuupenn/SpaGCN

Our tutorial: [run_SpaGCN](https://github.com/maiziezhoulab/BenchmarkST/blob/main/SpaGCN/run_SpaGCN/tutorial.ipynb)

#### DeepST

Github link: https://github.com/JiangBioLab/DeepST

Our tutorial: [run_DeepST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/DeepST/deepst/run_DeepST.ipynb)

#### BayesSpace

Github link: https://github.com/edward130603/BayesSpace

Our tutorial: [run_BayesSpace]()

#### BASS

Github link: https://github.com/zhengli09/BASS

Our tutorial: [run_BASS]()

#### PRECAST

Github link: https://github.com/cran/PRECAST

Our tutorial: [run_PRECAST]()

### original Github repository links and our script for Alignment task:

#### STAlign

Github link: https://github.com/JEFworks-Lab/STalign

Our tutorial: [run_STAlign]()

#### GPSA

Github link: https://github.com/andrewcharlesjones/spatial-alignment

Our tutorial: [run_GPSA]()

#### PASTE

Github link: https://github.com/raphael-group/paste

Our tutorial: [run_PASTE]()

#### PASTE2

Github link: https://github.com/raphael-group/paste2

Our tutorial: [run_PASTE2]()

#### SPACEL

Github link: https://github.com/QuKunLab/SPACEL

Our tutorial: [run_SPACEL]()

### original Github repository links and our script for Integration task:

#### SPIRAL

Github link: https://github.com/guott15/SPIRAL

Our tutorial: [run_SPIRAL]()

#### STAligner

Github link: https://github.com/zhanglabtools/STAligner

Our tutorial: [run_STAligner]()

#### PRECAST

Github link: https://github.com/cran/PRECAST

Our tutorial: [run_PRECAST]()

#### BASS

Github link: https://github.com/zhengli09/BASS

Our tutorial: [run_BASS]()

#### DeepST

Github link: https://github.com/JiangBioLab/DeepST

Our tutorial: [run_DeepST]()

## Citation

-------
Yunfei Hu, etc. "Benchmarking clustering, alignment, and integration methods for spatial transcriptomics", Genome Biology, under revision