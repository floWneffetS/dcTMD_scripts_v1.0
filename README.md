# dcTMD 

Python scripts used for dissipation corrected targeted molecular dynamics analysis. If you use them, please cite:
Wolf, S., Stock, G. Targeted Molecular Dynamics Calculations of Free Energy Profiles Using a Nonequilibrium Friction Correction. J. Chem. Theory Comput. 2018, DOI: 10.1021/acs.jctc.8b00835.

* NEQGamma.py: 
Integrates a constraint force file via trapezoid rule, calculates the NEQ memory friction kernel and friction factors, and performs a friction correction. First column: reaction coordinate in nm calculated via t * vel. Second column: force integral, i.e. the work profile. Third column: friction factors. Fourth column: trapezoid integral (final value) of friction work along reaction coordinate. Fourth column: friction corrected work profile. ATTENTION: Use with python3 or higher!

* NEQJarzy.py:
Integrates a constraint force file via trapezoid rule, and performs a friction correction based on Jarzynskis fast growth estimator (see Hendrix, D. A. and Jarzynski, C. J. Chem. Phys. 114, 5974 (2001)). First column: reaction coordinate in nm calculated via t * vel. Second column: mean work profile. Third column: friction factors derived via differentiation of dissipation work. Fourth column: Jarzynski estimator free energy profile. ATTENTION: Use with python3 or higher!
