# WUDESIM
A special water quality simulation model for the dead-end branches of drinking water distribution systems. This module is a C++ application that incorporates EPANET programmer's toolkit and solves the advection-dispersion-reaction (ADR) equation for solute transport in the links of the dead-end branches. The model generates a separate report (.rpt) file.  

- For more information:
> Abokifa, A. A., Yang, Y. J., Lo, C. S., & Biswas, P. (2016). Water quality modeling in the dead end sections of drinking water distribution networks. Water research, 89, 107-117.


- The compiled executable can be found in the \bin folder. The executable can be run directly from the command line. It takes four command line arguments:
 
> WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp WUDESIM.rpt

- EPANET.inp  ... Name of the EPANET input file
- EPANET.rpt  ... Desired name for the EPANET report file
- WUDESIM.inp ... Name of the corresponding WUDESIM input file
- WUDESIM.rpt ... Desired name for the WUDESIM report file

- An example of the four files can be found in the \test folder
