#!/bin/bash
CONDA_INSTALL_URL=${CONDA_INSTALL_URL:-"https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"}

####### Install conda locally ############
rm -rf conda_local/conda
rm -f /tmp/Miniconda3-latest-Linux-x86_64.sh
wget -P /tmp \
    "${CONDA_INSTALL_URL}" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p conda_local/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

export PATH=conda_local/conda/bin:$PATH
conda_local/conda/bin/python3 -m pip install nvidia-pyindex
conda env create --name=deepfoldrna -f environment.yml

######## Create environement for SPOT-RNA-1D #########
source scripts/activate_conda_env.sh 
conda create -n venv python=3.6
conda activate venv
conda install tensorflow==1.15.0
