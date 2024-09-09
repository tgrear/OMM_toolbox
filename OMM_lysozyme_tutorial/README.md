# OpenMM lysozyme in water tutorial

### 1) Create and activate OpenMM conda environment. This may require executing 'module load anaconda3' first if operating on the UNCC HPC.
    conda create -n OMM_env numpy matplotlib git
    conda activate OMM_env
    conda install -c conda-forge -c schrodinger pymol-bundle pdbfixer

### 2) Install OpenMM in conda environment and run OpenMM test:
    conda install openmm
    python -m openmm.testInstallation

### 3) Clone OMM_toolbox repository:
    git clone https://github.com/tgrear/OMM_toolbox.git
    cd OMM_toolbox/OMM_lysozyme_tutorial

### 4a) Run lysozyme in water simulation on UNCC HPC:
    sbatch runOMM.slurm

### 4b) If not on UNCC HPC, run lysozyme in water simulation locally:
    python OMM_lysozyme_tutorial.py
    
### The output is routed to the directory: OMM_output_DATE_TIME

---

### Additional OpenMM tutorials can be found [here](http://docs.openmm.org/latest/userguide/library/03_tutorials.html).

---
