# OMM_tools lysozyme in water tutorial

1) Create and activate OpenMM conda environment:\
    a) Execute: 'conda create -n OMM_env numpy matplotlib git'.\
    b) Execute: 'conda activate OMM_env'.\
    c) Execute: 'conda config --append channels conda-forge'.\
    d) Execute: 'conda install -c schrodinger pymol-bundle'.\

2) Install OpenMM in conda environment:\
    a) Execute: 'conda install openmm'.\
    b) Verify install: 'python -m openmm.testInstallation'.\

3) Clone OMM_tools lysozyme tutorial directory:\
    a) Execute: 'git clone https://github.com/tgrear/OMM_tools.git'.\
    b) Execute: 'cd OMM_tools/OMM_lysozyme_tutorial'.\

4) Run lysozyme in water simulation on UNCC HPC:\
    a) Execute: 'sbatch runOMM.slurm' if using UNCC HPC.\
    b) Otherwise, execute: 'python OMM_lysozyme_tutorial.py'.\
    c) The output is routed to the dir: 'OMM_output_DATE_TIME'.\

---
