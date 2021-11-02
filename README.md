# pmf_smd
make a potential of mean force (PMF) from several Steered Molecular Dynamics Simulation (SMD)

This soft computes the potential of mean force (PMF)  from several steered molecular dynamics (SMD) simulations. The PMF is computed through the second order cumulant method. Reference of this method could be found in:
"Free energy calculation from steered molecular dynamics simulations
using Jarzynski's equality". S. Park and F. Khalili-Araghi. 
J Chem Phys 2003 (119) 3559-66

Usage:
% pmf_smd.pl ^ basenamefile ^ number
 
For example, if you made 20 smd simulations and have files named: smd.txt.1 ; smd.txt.2 ... smd.txt.20. Then:

% pmf_smd.pl smd.txt 20

Several information will be asked to the user:
 - The starting value of the X-axis (reaction coordinate)
 - Ending value of the X-axis (reaction coordinate)
 - The number of PMF points

Beware: starting value < ending value
#
# Several files will be created:
# + PMF_curve.dat
# + PMF_gnuplot.gnu
# + PMF_number.png
# 
# PMF_curve.dat is the most important file. This one contains 5 columns
# which represent, respectively, the reaction coordinate, the free energy
# with the first-order cumulant method, the error bar of this last value,
# the free energy with the second-order cumulant method, the error bar of
# this last value and the last column is the Herman error (J Phys Chem 
# 1991 (95) 9029-34. PMF_gnuplot.gnu is a gnuplot script file which is used
# to provide the PMF_number.png.
#
# Enjoy, Florent Barbault
