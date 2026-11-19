
# RenormHNet
RenormHNet (**renorm**alized **H**amiltonian neural **net**work) is the code released with the paper "Machine-Learning Spin-Fermion Models for Strongly Disordered Magnets: Sb-Doped MnTe". 

This code works as a patch of [DeepH-pack](https://github.com/mzjb/DeepH-pack.git), 
and preserves the original workflow of DeepH. The installation and use methods are almost the same. 

The code newly includes the embedding of spin Hamiltonian with Heisenberg-type exchange 
couplings obtained from [Jx](https://kaist-elst.github.io/Jx.jl/) (DOI: https://doi.org/10.1016/j.cpc.2019.106927). 
The currently supported DFT package is [OpenMX](http://www.openmx-square.org/). 


# Contents
1. [Installation](#installation)
1. [Usage](#usage)
1. [Contributors](#contributors)
1. [howtocite](#how-to-cite)


# Installation
RenormHNet requires Python and Julia. 

## Python
Prepare a Python 3.9 interpreter. 

To set a conda environment on Linux:
```bash 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
If you already have conda installed: 
```bash
conda create -n renormhnet python=3.9.21
conda activate renormhnet

conda install numpy scipy pip -y

conda install pytorch==1.9.1 ${pytorch_config}
# ${pytorch_config} is replaced by your configuration. For example,
#conda install pytorch==1.9.1 cudatoolkit=11.1 -c pytorch && python -m pip install --no-cache-dir torchvision==0.10.1 torchaudio==0.9.1 && python -m pip install --no-cache-dir --only-binary=:all: torch-scatter==2.0.9 torch-sparse==0.6.12 -f https://data.pyg.org/whl/torch-1.9.1+cu111.html

python -m pip install torch-geometric==1.7.2
conda install pymatgen -c conda-forge
python -m pip install opt_einsum_fx==0.1.4
python -m pip install e3nn==0.3.5
python -m pip install h5py
python -m pip install tensorboard
python -m pip install pathos
python -m pip install psutil

conda install pymatgen -c conda-forge
```

For the pytorch configuration, se the [official PyTorch website](https://pytorch.org/get-started/previous-versions/) for how to set it. 

## Julia
The Julia versions 1.11.3, 1.12.6, and 1.13.4 have been tested. Higher versions are also expected to work.
To install Julia, visit https://julialang.org/.
```
wget https://julialang-s3.julialang.org/bin/linux/x64/1.13/julia-1.13.4-linux-x86_64.tar.gz
tar -zxvf julia-1.13.4-linux-x86_64.tar.gz
```
Add the following to your `~/.bashrc` or `~/.profile`:
```
export PATH="$PATH:${JULIAPATH}/bin"
```
${JULIAPATH} is the path you've just installed.

We recommend creating an independent environment for RenormHNet using ``julia --project=${PATH_TO_BE_INSTALLED}``.

Run ``julia`` in a terminal, then type `]` to enter the Pkg REPL. 

```julia
pkg> add Arpack HDF5 ArgParse JLD JSON IterativeSolvers DelimitedFiles StaticArrays LinearMaps DataFrames CSV
```

## Jx 
Please refer to the official [Jx installation guide](https://kaist-elst.github.io/Jx.jl/docs/quick-start-guide/).
```bash
git clone https://github.com/KAIST-ELST/Jx.jl
cd Jx.jl
julia install.jl
```

## DFT package
[OpenMX3.9](http://www.openmx-square.org/download.html) is currently the only supported DFT package.  
As in DeepH-pack, the overlap-only version of OpenMX is also required; install it from its [repository](https://github.com/mzjb/overlap-only-OpenMX).

## Install RenormHNet
```bash
git clone https://github.com/DHKiem/RenormHNet.git
cd RenormHNet
conda activate renormhnet
pip install .
```

# Usage

# Contributors

# How to cite