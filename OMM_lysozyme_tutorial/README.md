# OpenMM lysozyme in water tutorial

### 1) Create and activate OpenMM conda environment:
    conda create -n OMM_env numpy matplotlib git 
    conda activate OMM_env 
    conda config --append channels conda-forge 
    conda install -c schrodinger pymol-bundle 

### 2) Install OpenMM in conda environment:
    conda install openmm
    python -m openmm.testInstallation

### 3) Clone OMM_tools repository:
    git clone https://github.com/tgrear/OMM_tools.git
    cd OMM_tools/OMM_lysozyme_tutorial

### 4) Run lysozyme in water simulation on UNCC HPC:
    a) Execute: 'sbatch runOMM.slurm' if using UNCC HPC.\
    b) Otherwise, execute: 'python OMM_lysozyme_tutorial.py
    c) The output is routed to the dir: 'OMM_output_DATE_TIME

---
