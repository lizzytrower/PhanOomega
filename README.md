# PhanOomega

This suite of Matlab code is designed to reconstruct the carbonate chemistry of Phanerozoic seawater using ooid size data. The main codes and the order in which to run them are:

1. PhanOomega_ooid_size_data_input.m: this code loads data from a spreadsheet (so could be used on a filtered or expanded dataset), organizes the data, and generates some plots that are useful for determining if the data read in correctly. It generates a simplified output file that can be used for the main model. This code also applies a size correction (for size data estimated from thin sections) using the ooidsizecorrection.m function file.
2. PhanOomega_model.m: this is the main model code that calculates Omega, then DIC, Alk, and pH based on the ooid size data. This code uses a number of custom functions, including aragoniteinterp.m, calciteinterp_BW.m, and calciteinterp_L.m (which interpolate precipitation kinetics based on temperature), Oomegasolver.m (which estimates Omega based on ooid diameter and requires the additional file susp_abrasion_calculations_abrcalc.m to run), and phreeqc_oomega_aragonite.m and phreeqc_oomega_calcite.m (which run a COM server version of PHREEQC - must be installed separately - to calculate carbonate speciation.
3. PhanOomega_mainfigures.m: this code generates the main figures in the accompanying manuscript
4. PhanOomega_suppfigures: this code generates supplementary figures S2 and S5 in the accompanying manuscriot

External code needed:
For plotting: these codes use the magma and viridis perceptually uniform Matlab colormaps, which can be downloaded here: https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps
For PHREEQC: the main code requires the PHREEQC COM server. Downloads and documentation can be found here: https://www.usgs.gov/software/phreeqc-version-3
