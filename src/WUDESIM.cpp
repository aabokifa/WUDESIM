/*
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
*/

#include <iostream>
#include <ctime>

# include "WUDESIM.h"
# include "CLASSES.h"
# include "WUDESIM_CORE.h"
# include "WRITING_FUN.h"

// Define network
Network net;

// initialize error
int    error = 0;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Run complete simulation

int DE_RUN_FULL_SIM(char* EPANET_INP, char* EPANET_RPT, string WUDESIM_INP, string WUDESIM_RPT) {

	error = DE_OPEN_EPANET_PROJ(EPANET_INP); if (error) { return 1; }

	error = DE_FIND_DEADENDS(); if (error) { return 1; }

	error = DE_RUN_EPANET_SIM(EPANET_RPT); if (error) { return 1; }

	error = DE_CALC_DEADEND_PROPERTIES(); if (error) { return 1; }

	error = DE_OPEN_WUDESIM_PROJ(WUDESIM_INP); if (error) { return 1; }

	error = DE_GENERATE_STOC_DEMAND(); if (error) { return 1; }

	error = DE_RUN_WUDESIM_SIM(WUDESIM_RPT); if (error) { return 1; }

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Open EPANET Project
int DE_OPEN_EPANET_PROJ(char* EPANET_INP) {

	// Get file names	
	net.EPANET_INP = EPANET_INP;

	// Process EPANET input file
	cout << "\no	Processing EPANET input file  ..." << endl;
	error = OP_EPANET_INP(&net);
	if (error) { cout << "o	Processing EPANET input file was not successful!" << endl; return 1; }
	else { cout << "o	Processing EPANET input file was successful!" << endl; }

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Find dead-end branches in the network
int DE_FIND_DEADENDS() {
	
	// Find dead end branches in the network
	cout << "\no	Finding dead-end branches ..." << endl;
	int N_DEbranches = 0;
	N_DEbranches = FIND_DE_BRANCHES(&net);
	if (N_DEbranches == 0) { cout << "o	Found no dead-end branches in the network!" << endl; return 1; }
	else { 
		cout << "o	Found " << N_DEbranches << " dead-end branches in the network -> See DE_Pipe_ID.out" << endl; 
		write_DE_ids(&net); 
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Run EPANET simulation
int DE_RUN_EPANET_SIM(char* EPANET_RPT) {
	
	// Get file names
	net.EPANET_RPT = EPANET_RPT;

	// Run EPANET simulation and extract results
	cout << "\no	Running EPANET simulation..." << endl;
	error = RUN_EPANET_SIM(&net);
	if (error) { cout << "o	Running EPANET was not successful!" << endl; return 1; }
	else { cout << "o	Running EPANET was successful -> See EPANET report file!" << endl; }

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate the properties of dead-end branches
int DE_CALC_DEADEND_PROPERTIES(){

	//   Calculate dead end pipes properties
	cout << "\no	Evaluating properties of dead-end pipes ..." << endl;
	error = CALC_DE_PROPERTIES(&net);
	if (error) { cout << "o	Evaluating properties of dead-end pipes was not successful!" << endl; return 1; }
	else { 
		cout << "o	Evaluating properties of dead-end pipes was successful -> See DE_Properties.out" << endl; 
		write_DE_Properties(&net);
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Open WUDESIM Project
int DE_OPEN_WUDESIM_PROJ(string WUDESIM_INP) {

	// Get file names	
	net.WUDESIM_INP = WUDESIM_INP;

	// Process WUDESIM input file
	cout << "\no	Processing WUDESIM input file ..." << endl;
	error = OP_WUDESIM_INP(&net);
	if (error) { cout << "o	Processing WUDESIM input file was not successful!" << endl; return 1; }
	else { cout << "o	Processing WUDESIM input file was successful!" << endl; }

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Generate stochastic demands for deadend pipes
int DE_GENERATE_STOC_DEMAND() {
	
	//  Generate stochastic demands
	if (net.DE_options.Stoc_dem_fl) {
		cout << "\no	Generating stochastic flows for dead-end pipes ..." << endl;
		error = GEN_STOC_DEM(&net);
		if (error) { cout << "o	Stochastic flow generation was not successful!" << endl; return 1; }
		else { 
			cout << "o	Stochastic flow generation was successful -> See DE_Stochastic_Flow.out" << endl; 
			write_stoc_dems(&net);
		}	
	}

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Run WUDESIM simulation
int DE_RUN_WUDESIM_SIM(string WUDESIM_RPT) {

	// Get file name
	net.WUDESIM_RPT = WUDESIM_RPT;

	cout << "\no	Starting WUDESIM water quality simulations!" << endl;
	error = RUN_WUDESIM_SIM(&net);
	if (error) { cout << "o	WUDESIM simulations were not successful!" << endl; return 1; }
	else { 
		cout << "o	WUDESIM simulations were successful -> See WUDESIM report file!" << endl; 
		write_WUDESIM_rpt(&net);
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////