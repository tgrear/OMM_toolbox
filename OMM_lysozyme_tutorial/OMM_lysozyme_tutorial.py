###########################################################################################
#                                   |                                                    ##
# File:    OMM_lysozyme_tutorial.py |             BioMolecular Physics Group             ##
# Author:  Tyler J Grear            |            University of North Carolina            ##
# Created: 07-18-2024               |                    at Charlotte                    ##
#                                   |                                                    ##
###########################################################################################
#>=================================<| Import Libraries |>================================<#

import os                        # Import working operating system module.
import sys                       # System-specific parameters and functions.
import time                      # Import module to benchmark run times.
import datetime                  # Import time module for time-related functions.
import subprocess                # Import subprocess module.
import numpy as np               # Import Numerical Python library.
import typlot as tlt             # Import typlot package for custom visualizations.
from openmm import *             # Import all packages from openmm library.
from openmm.app import *         # Import all packages from openmm.app library.
from openmm.unit import *        # Import all packages from openmm.unit library.
import matplotlib.pyplot as plt  # Import standard plotting library.

#>=======================================================================================<#
#>====================================<| Parameters |>===================================<#

pName = "1AKI"        # Set name of protein that corresponds to pdb file (i.e., 1AKI.pdb)
FF = "charmm36"       # Charmm36 ("charmm36") force field
H2O = "water"         # Charmm36 ("water") water model
temp = 300            # Set temperature (in Kelvin) for simulation
n_md = 10000          # Set the number of molecular dynamics simulation steps
n_rep = 100           # Set report/logging rate, n_rep/step
padding = 20          # Set padding for simulation box.
oName = "OMM_output"  # Set name of output directory
HPC = False           # Set to True if using CUDA for GPUs, False will bypass this feature

#>=======================================================================================<#
#>================================<| Simulation Config |>================================<#

DT = datetime.datetime.now().strftime("%m-%d-%Y_%H-%M-%S")  # Get run date-time (DT)
out_dir = oName+"_"+DT+os.sep         # Set output directory name
if os.path.exists(out_dir) == False:  # If out directory !exist
   os.mkdir(out_dir)                  # Create output directory

if HPC == True:
    platform = Platform.getPlatformByName("CUDA")

print(""); print("Executing OpenMM..."); print("System: "+pName); print("")
start = time.time(); start0 = time.time()
pdb = PDBFile("data"+os.sep+pName+".pdb")
forcefield = ForceField(FF+".xml",FF+os.sep+H2O+".xml")
modeller = Modeller(pdb.topology,pdb.positions)
print("Removing crystallized water molecules...")
modeller.deleteWater()
print("Adding protonation states and missing hydrogen atoms...")
modeller.addHydrogens(forcefield)

positions = np.array(modeller.positions.value_in_unit(angstrom))
box_center = (np.max(positions,axis=0) - np.min(positions,axis=0) + padding)/2
box_vectors = [vec3.Vec3(box_center[0]*2,0.0,0.0)*angstrom,
               vec3.Vec3(0.0,box_center[1]*2,0.0)*angstrom,
               vec3.Vec3(0.0,0.0,box_center[2]*2)*angstrom]

means = np.mean(np.array(modeller.positions.value_in_unit(angstrom)),axis=0)
means -= box_center; c = vec3.Vec3(*means)
modeller.positions = [p.value_in_unit(angstrom) - c for p in modeller.positions]*angstrom

modeller.addSolvent(forcefield,padding=1.0*nanometer)
system = forcefield.createSystem(modeller.topology,
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1.0*nanometer,
                                 constraints=HBonds)
integrator = LangevinMiddleIntegrator(temp*kelvin,1/picosecond,0.004*picoseconds)
simulation = Simulation(modeller.topology,system,integrator)
simulation.context.setPositions(modeller.positions)
print("Simulation setup complete: "+str(time.time() - start)+" seconds."); print("")

#>=======================================================================================<#
#>===============================<| Energy Minimization |>===============================<#

print("Minimizing energy...")
start = time.time()
simulation.minimizeEnergy()
print("Energy minimization complete: "+str(time.time() - start)+" seconds."); print("")

#>=======================================================================================<#
#>===================================<| Sim Logging |>===================================<#

simulation.reporters.append(PDBReporter(out_dir+"sim_"+pName+".pdb",n_rep))
simulation.reporters.append(StateDataReporter(sys.stdout,n_rep,
                                              step=True,
                                              potentialEnergy=True,
                                              kineticEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              volume=True,
                                              density=True,
                                              systemMass=None))
simulation.reporters.append(StateDataReporter(out_dir+"md_log.txt",
                                              int(n_rep/10),
                                              step=True,
                                              potentialEnergy=True,
                                              kineticEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              volume=True,
                                              density=True,
                                              systemMass=None))

#>=======================================================================================<#
#>=====================================<| NVT Eq. |>=====================================<#

print("Running NVT equillibration..."); start = time.time()
simulation.step(n_md)
print("NVT equillibration complete: "+str(time.time() - start)+" seconds."); print("")

#>=======================================================================================<#
#>=====================================<| NPT Eq. |>=====================================<#

system.addForce(MonteCarloBarostat(1*bar,temp*kelvin))
simulation.context.reinitialize(preserveState=True)
print("Running NPT equillibration..."); start = time.time()
simulation.step(n_md)
print("NPT equillibration complete: "+str(time.time() - start)+" seconds."); print(" ")

#>=======================================================================================<#
#>=================================<| Save Simulation |>=================================<#

s = simulation.context.getState(getPositions=True,
                                getVelocities=True,
                                getEnergy=True,
                                getForces=True)
with open(out_dir+"state_"+pName+".xml","w") as f: f.write(XmlSerializer.serialize(s))
print(pName+"_sim.pdb saved to: "+out_dir)
print(pName+"_state.xml saved to: "+out_dir); print("")

#>=======================================================================================<#
#>==================================<| Visualization |>==================================<#

tlt.simPlot(pName,out_dir)
pout = out_dir+"pymol"+os.sep
tlt.genMv(pName,out_dir+"sim_"+pName+".pdb",pout)

#>=======================================================================================<#

print("OpenMM finished: "+str(time.time() - start0)+" seconds.")

#>=======================================================================================<#