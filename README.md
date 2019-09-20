# WUDESIM

A public-domain software package for simulating water quality in the dead-end branches of drinking water distribution networks. The C/C++ source code incorporates EPANET programmers' toolkit, and solves the advection-dispersion-reaction (ADR) equation for solute transport in the pipes of the dead-end branches. The software generates a separate report (.RPT) file together with a series of output (.OUT) files. 

- The compiled Windows executable and DLL toolkit can be found in the \bin folder. The executable can be run directly from the command line and takes four command line arguments:
 
> WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp WUDESIM.rpt

- EPANET.INP  ... Name of the EPANET input file
- EPANET.RPT  ... Desired name for the EPANET report file
- WUDESIM.INP ... Name of the corresponding WUDESIM input file
- WUDESIM.RPT ... Desired name for the WUDESIM report file

Examples on how the users can employ the executable as well as the different toolkit functions to run water quality analysis and obtain simulation results are provided in the \test folder
