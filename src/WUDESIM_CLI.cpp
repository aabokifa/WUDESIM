/* 
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
*/

#include <iostream>

#include "WUDESIM.h"

using namespace std;

int main(int argc, char* argv[])
{
	// Check number of input CL arguments is correct
	if (argc != 5) {
		cout << "Command line input should be: WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp WUDESIM.rpt" << endl;
		return 1;
	}

	char * EPANET_INP  = argv[1];
	char * EPANET_RPT  = argv[2];
	char * WUDESIM_INP = argv[3];
	char * WUDESIM_RPT = argv[4];

	// Run a complete analysis
	int error = DE_RUN_FULL_SIM(EPANET_INP, EPANET_RPT, WUDESIM_INP, WUDESIM_RPT);

	if (error) { return 1; }

	return 0;
}