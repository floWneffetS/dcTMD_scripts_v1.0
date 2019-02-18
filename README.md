# dcTMD 

Python scripts used for dissipation corrected targeted molecular dynamics analysis for usage with "\*pullf.xvg" files from Gromacs v2016 or higher (Application of usage with Gromacs v5.X and lower only needs adjustment in the "skiprows" argument in the data read in step). If you use them, please cite:
Wolf, S., Stock, G. Targeted Molecular Dynamics Calculations of Free Energy Profiles Using a Nonequilibrium Friction Correction. J. Chem. Theory Comput. 2018, 14(12), 6175-6182, DOI: 10.1021/acs.jctc.8b00835.

* NEQGamma.py: 
Integrates a constraint force file via trapezoid rule, calculates the NEQ memory friction kernel and friction factors, and performs a friction correction. Use "-h" for more information on usage.
ATTENTION: Use with python3 or higher!

* NEQJarzy.py:
Integrates a constraint force file via trapezoid rule, and performs a friction correction based on Jarzynskis fast growth estimator (see Hendrix, D. A. and Jarzynski, C. J. Chem. Phys. 114, 5974 (2001)). Use "-h" for more information on usage. 
ATTENTION: Use with python3 or higher!
