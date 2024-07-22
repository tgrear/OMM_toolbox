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

### 4a) Run lysozyme in water simulation on UNCC HPC:
    sbatch runOMM.slurm

### 4b) Run lysozyme in water simulation locally:
    python OMM_lysozyme_tutorial.py
    
### The output is routed to the directory: 'OMM_output_DATE_TIME

---
