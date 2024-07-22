###########################################################################################

import os                        # Import working operating system module.
import time                      # Import module to benchmark run times.
import pymol                     # Import pymol package for MDS visualization.
import datetime                  # Import time module for time-related functions.
import subprocess                # Import subprocess module.
import numpy as np               # Import Numerical Python library.
import matplotlib.pyplot as plt  # Import standard plotting library.

#>=======================================================================================<#
#>-------------------------------------<| simPlot |>-------------------------------------<#
#                                                                                         #
# INPUT:    1) pName   | Name of protein (i.e., 1AKI).                                    #
#           2) out_dir | Path of output directory to generate simulation figures.         #
#                                                                                         #
#>---------------------------------------------------------------------------------------<#

def simPlot(pName,out_dir):
    print("Generating figures..."); start = time.time()
    data = np.loadtxt(out_dir+os.sep+"md_log.txt", delimiter=',')
    step = data[:,0]; PE = data[:,1]; KE = data[:,2]; Etot = data[:,3]
    T = data[:,4]; V = data[:,5]; D = data[:,6]
    #>====================================================================<#
    plt.plot(step,PE,color="black",zorder=6,linewidth=1)
    plt.axhline(y=np.mean(PE),color="red",linestyle=":",linewidth=2,zorder=5) 
    plt.xlabel("Step"); plt.ylabel("Potential energy (kJ/mol)")
    plt.title(pName+r": PE vs. step  |  $\mu_{PE}$ = "+
              str(round(np.mean(PE),4))+" kJ/mol")
    plt.legend([r"PE",r"$\mu_{PE}$"],loc='best',edgecolor="k",ncol=2)
    plt.grid(); plt.show(); plt.tight_layout(pad=0.165)
    plt.savefig(out_dir+"PE_vs_step.png",bbox_inches='tight',dpi=300)
    plt.close('all')
    #>====================================================================<#
    plt.plot(step,KE,color="black",zorder=6,linewidth=1)
    plt.axhline(y=np.mean(KE),color="red",linestyle=":",linewidth=2,zorder=5) 
    plt.xlabel("Step"); plt.ylabel("Kinetic energy (kJ/mol)")
    plt.title(pName+r": KE vs. step  |  $\mu_{KE}$ = "+
              str(round(np.mean(KE),4))+" kJ/mol")
    plt.legend([r"KE",r"$\mu_{KE}$"],loc='best',edgecolor="k",ncol=2)
    plt.grid(); plt.show(); plt.tight_layout(pad=0.165)
    plt.savefig(out_dir+"KE_vs_step.png",bbox_inches='tight',dpi=300)
    plt.close('all')
    #>====================================================================<#
    plt.plot(step,Etot,color="black",zorder=6,linewidth=1)
    plt.axhline(y=np.mean(Etot),color="red",linestyle=":",linewidth=2,zorder=5) 
    plt.xlabel("Step"); plt.ylabel("Total energy (kJ/mol)")
    plt.title(pName+r": $\Sigma$E vs. step  |  $\mu_{\Sigma E}$ = "+
              str(round(np.mean(Etot),4))+" kJ/mol")
    plt.legend([r"$\Sigma$E",r"$\mu_{\Sigma E}$"],loc='best',edgecolor="k",ncol=2)
    plt.grid(); plt.show(); plt.tight_layout(pad=0.165)
    plt.savefig(out_dir+"Etot_vs_step.png",bbox_inches='tight',dpi=300)
    plt.close('all')
    #>====================================================================<#
    plt.plot(step,T,color="black",zorder=6,linewidth=1)
    plt.axhline(y=np.mean(T),color="red",linestyle=":",linewidth=2,zorder=5) 
    plt.xlabel("Step"); plt.ylabel(r"Temperature (K$\degree$)")
    plt.title(pName+r": Temperature vs. step  |  $\mu_{T}$ = "+
              str(round(np.mean(T),4))+r" K$\degree$")
    plt.legend([r"T",r"$\mu_{T}$"],loc='best',edgecolor="k",ncol=2)
    plt.grid(); plt.show(); plt.tight_layout(pad=0.165)
    plt.savefig(out_dir+"T_vs_step.png",bbox_inches='tight',dpi=300)
    plt.close('all')
    #>====================================================================<#
    plt.plot(step,V,color="black",zorder=6,linewidth=1)
    plt.axhline(y=np.mean(V),color="red",linestyle=":",linewidth=2,zorder=5) 
    plt.xlabel("Step"); plt.ylabel(r"Box volume (nm$^{3}$)")
    plt.title(pName+r": Box volume vs. step  |  $\mu_{V}$ = "+
              str(round(np.mean(V),4))+r" nm$^{3}$")
    plt.legend([r"V",r"$\mu_{V}$"],loc='best',edgecolor="k",ncol=2)
    plt.grid(); plt.show(); plt.tight_layout(pad=0.165)
    plt.savefig(out_dir+"V_vs_step.png",bbox_inches='tight',dpi=300)
    plt.close('all')
    #>====================================================================<#
    plt.plot(step,D,color="black",zorder=6,linewidth=1)
    plt.axhline(y=np.mean(D),color="red",linestyle=":",linewidth=2,zorder=5) 
    plt.xlabel("Step"); plt.ylabel(r"Density (g/mL)")
    plt.title(pName+r": Density vs. step  |  $\mu_{\rho}$ = "+
              str(round(np.mean(D),4))+r" g/mL")
    plt.legend([r"$\rho$",r"$\mu_{\rho}$"],loc='best',edgecolor="k",ncol=2)
    plt.grid(); plt.show(); plt.tight_layout(pad=0.165)
    plt.savefig(out_dir+"D_vs_step.png",bbox_inches='tight',dpi=300)
    plt.close('all')
    #>====================================================================<#
    print("Simulation plots generated: "+str(time.time() - start)+" seconds."); print(" ") 

#>=======================================================================================<#
#>--------------------------------------<| genMv |>--------------------------------------<#
#                                                                                         #
# INPUT:    1) pName | Name of protein (i.e., 1AKI).                                      #
#           2) fname | Name of output pdb from simulation (i.e., 1AKI_sim.pdb).           #
#           3) pout  | Path of output directory to generate pymol sessions.               #
#                                                                                         #
#>---------------------------------------------------------------------------------------<#

def genMv(pName,fName,pout):
    print("Launching pymol..."); start = time.time()
    if os.path.exists(pout) == False:  # If pymol output directory does not exist
       os.mkdir(pout)                  # Create pymol output directory
    subprocess.run(["pymol","-qcQ",fName,  # Launch cli pymol on quiet mode with no splash
        "-d","zoom complete=1",            # Zoom to inculde entire PBC box
        "-d","remove resn CL",             # Remove chloride ions
        "-d", "preset.ball_and_stick(selection='all')",        # Set preset configuration
        "-d", "save "+pout+pName+"_ball_and_stick_wSolv.pse",  # Save animation
        #>===============================================================================<#
        "-d", "remove solvent",                          # Remove solvent
        "-d", "show cell",                               # Show PBC simulation box
        "-d", "save "+pout+pName+"_ball_and_stick.pse",  # Save animation
        #>===============================================================================<#
        "-d", "set bg_rgb,[1,1,1]",           # Set background color to white
        "-d", "color black, (name C*)",       # Set color of carbons to black
        "-d", "show surface, All",            # Show surface
        "-d", "set transparency, 0.5",        # Set surface transparency to 0.5
        "-d", "save "+pout+pName+"_TG1.pse",  # Save animation
        #>===============================================================================<#
        "-d", "set transparency, 0",          # Set surface transparency to 0
        "-d", "save "+pout+pName+"_TG2.pse",  # Save animation
        #>===============================================================================<#
        "-d", "set transparency, 1",             # Set surface transparency to 1
        "-d", "hide sticks",                     # Hide sticks representation
        "-d", "hide spheres",                    # Hide spheres representation
        "-d", "show cartoon, All",               # Show cartoon for all
        "-d", "color firebrick, ss h",           # Set helix color to firebrick
        "-d", "color deepblue, ss s",            # Set sheet color to deepblue
        "-d", "color grey70, ss l+''",           # Set loops color to grey70
        "-d", "set cartoon_discrete_colors, 1",  # Hide sticks representation
        "-d", "save "+pout+pName+"_TG3.pse",     # Save animation
        #>===============================================================================<#
        "-d", "preset.technical(selection='all')",  # Set preset configuration
        "-d", "show cell",                          # Show PBC simulation box
        "-d", "set bg_rgb,[0,0,0]",                 # Set background color to black
        "-d", "save "+pout+pName+"_technical.pse",  # Show PBC simulation box
        #>===============================================================================<#
        "-d", "quit"])  # Quit pymol cli
    print("Pymol animations generated: "+str(time.time() - start)+" seconds."); print(" ")

###########################################################################################