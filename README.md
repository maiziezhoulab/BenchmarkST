# BenchmarkST

We present a comprehensive evaluation of the performance of state-of-the-art methods for cell/spot clustering on spatial transcriptomic datasets. Our goal is to provide the research field with a fair evaluation and comparison of all available methods to facilitate the users’ choice as well as shed more light on further developments to focus on the challenging aspects of ST clustering task.

## Table of Content 
- [Installation](#installation)
- [Repository description](#Repository-description)
- [General Usage](#General-Usage)
  -  [ADEPT](#ADEPT)
  -  [GraphST](#GraphST)
  -  [conST](#conST)
  -  [WGS mode Small Indel detection](#wgs-mode-small-Indel-detection)
- [Citation](#Citation)


## Installation

We strongly suggest to follow the installation instructions given by each method. We will also provide the .yaml file of our debugging environment in each subfolder.

## Repository description
We provide all the jupyter notebooks in each subfolder to run and evaluate all algorithms for all datasets with ground truth annotation, and to reproduce the results introduced in our paper.


1. '${algorithm_name}' folders contain mostly the original repo's content and code, but we also adapt a bit so that it could be run on the new datasets that have not been tested in their published paper.
2. 'analysis' folder constains script for reproducing the figures in our paper.
3. -------
4. -------

## General Usage

Please follow the readme or jupyter notebook instructions in each subfolder for further details.

### original Github repository links and our script running tutorials:

#### ADEPT

Github link: https://github.com/maiziezhoulab/AquilaDeepFilter

Our tutorial: [run_ADEPT](https://github.com/maiziezhoulab/BenchmarkST/blob/main/ADEPT/run_ADEPT/run_ADEPT.ipynb)


#### conST

Github link: https://github.com/ys-zong/conST

Our tutorial: [run_conST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/conST/run_conST.ipynb)


#### GraphST

Github link: https://github.com/JinmiaoChenLab/GraphST

installation: https://deepst-tutorials.readthedocs.io/en/latest/Installation.html

Our tutorial: [run_GraphST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/GraphST/run_GraphST/tutorial.ipynb)

#### DeepST

Github link: https://github.com/JiangBioLab/DeepST

Our tutorial: [run_DeepST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/DeepST/deepst/run_DeepST.ipynb)

#### conGI

Github link: https://github.com/biomed-AI/ConGI

Our tutorial: [run_conGI](https://github.com/maiziezhoulab/BenchmarkST/blob/main/ConGI/run_ConGI/run_conGI.ipynb)

#### STAGATE

Github link: https://github.com/QIFEIDKN/STAGATE_pyG

installation: https://stagate.readthedocs.io/en/latest/Installation_pyG.html

Our tutorial: [run_STAGATE](https://github.com/maiziezhoulab/BenchmarkST/blob/main/STAGATE_pyG/run_STAGATE/run_STAGATE.ipynb)

#### CCST

Github link: https://github.com/xiaoyeye/CCST

Our tutorial: [run_CCST](https://github.com/maiziezhoulab/BenchmarkST/blob/main/CCST/run_CCST.ipynb)

#### SpaGCN

Github link: https://github.com/jianhuupenn/SpaGCN

Our tutorial: [run_SpaGCN](https://github.com/maiziezhoulab/BenchmarkST/blob/main/SpaGCN/run_SpaGCN/tutorial.ipynb)

#### SEDR

Github link: https://github.com/JinmiaoChenLab/SEDR

Our tutorial: [run_SEDR](https://github.com/maiziezhoulab/BenchmarkST/blob/main/SEDR/run_SEDR.ipynb)

## Citation

-------