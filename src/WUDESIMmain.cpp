/* 
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: A software for water quality simulation in the dead-end sections of drinking water distribution systems
*/

#include <iostream>
#include <ctime>

#include "Classes.h"
#include "WUDESIMmain.h"
#include "Utilities.h"
#include "epanet2.h"

using namespace std;

int main(int argc, char* argv[])
{
	int start_s = clock();
	int error = 0, N_DEbranches = 0, N_xjuncts = 0;

	//Check number of input CL arguments is correct
	if (argc != 5) {
		cout << "Command line input should be: WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp WUDESIM.rpt" << endl;
		return 1;
	}
	else {
		cout << "****************	WUDESIM Started!	****************" << endl;
	}
		
	char *INPfileName, *RPTfileName, *WUDESIMinp, *WUDESIMrpt;

	INPfileName = argv[1];
	RPTfileName = argv[2];
	WUDESIMinp  = argv[3];
	WUDESIMrpt  = argv[4];

	// Define network
	Network net;

	// Process EPANET input file
	error = OpenEPANETinp(INPfileName, &net);
	if (!error) {
		cout << "o	Processing EPANET input file was successful" << endl;
	}
	else {
		cout << "o	Processing EPANET input file was not successful" << endl;
		return 1;
	}


	// Find dead end branches in the network
	N_DEbranches = DEFIND(&net);
	if (N_DEbranches) {
		cout << "o	Found " << N_DEbranches << " dead-end branches in the network --> See DEpipeID.txt" << endl;
	}
	else {
		cout << "o	Found no dead-end branches in the network" << endl;
		return 1;
	}

	/*
	// Find cross junctions in the network
	N_xjuncts = XJFIND(&net);
	cout << "o	Found " << N_xjuncts << " cross junctions in the network --> See XjuncIDs.txt" << endl;
	*/

	// Run EPANET simulation and extract results
	error = OpenEPANETrpt(INPfileName, RPTfileName, &net);
	if (!error) {
		cout << "o	Running EPANET was successful" <<endl;
		cout << "o	Starting WUDESIM simulations" << endl;
		cout << "*****************************************************************************" << endl;
	}
	else {
		cout << "o	Running EPANET was not successful" << endl;
		return 1;
	}

	//  Run Water Quality simulations for dead-end branches
	error = WQSIM(WUDESIMinp, WUDESIMrpt, &net);
	if (!error) {
		cout << "o	WUDESIM simulations finished successfuly" << endl;
		cout << "****************	Exiting WUDESIM!	****************" << endl;

	}
	else {
		cout << "o	WUDESIM simulations were not successful" << endl;
		return 1;
	}

	int stop_s = clock();
	cout << "execution time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000 << " ms" << endl;
	return 0;
}