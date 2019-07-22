/* 
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
*/

#include <iostream>
#include <ctime>

#include "WUDESIM.h"

using namespace std;

int main(int argc, char* argv[])
{
	// Start clock
	double start_s = clock();

	// Check number of input CL arguments is correct
	if (argc != 5) {
		cout << "Command line input should be: WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp WUDESIM.rpt" << endl;
		return 1;
	}
	else {
		cout << "****************	WUDESIM Started!	****************" << endl;
	}

	char * EPANET_INP  = argv[1];
	char * EPANET_RPT  = argv[2];
	string WUDESIM_INP = argv[3];
	string WUDESIM_RPT = argv[4];

	// Run a complete analysis
	int error = DE_RUN_FULL_SIM(EPANET_INP, EPANET_RPT, WUDESIM_INP, WUDESIM_RPT);

	if (error) { cout << "WUDESIM did not finish successfuly!" << endl; return 1; }

	else {
		double stop_s = clock();
		cout << "\n*****************************************************************************" << endl;
		cout << "WUDESIM finished successfuly!" << endl;
		cout << "WUDESIM execution time: " << stop_s - start_s / double(CLOCKS_PER_SEC) * 1000 << " ms" << endl;
	}
	return 0;
}