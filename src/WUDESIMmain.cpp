/* 
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: A software for water quality simulation in the dead-end sections of drinking water distribution systems
*/

#include <iostream> 

#include "Classes.h"
#include "WUDESIMmain.h"
#include "Utilities.h"

using namespace std;

int main(int argc, char *argv[])
{
	int error = 0, N_DEbranches = 0, N_xjuncts = 0;

	//Check number of input CL arguments is correct
	if (argc != 4) {
		cout << "Command line input should be: WUDESIM.exe EPANET.inp EPANET.rpt WUDESIM.inp" << endl;
		return 1;
	}


	char *INPfileName, *RPTfileName, *WUDESIMinp;

	INPfileName = argv[1];
	RPTfileName = argv[2];
	WUDESIMinp  = argv[3];

	// Define network
	Network net;

	// Process EPANET input file
	error = OpenEPANETinp(INPfileName, &net);
	if (!error) {
		cout << "Processing EPANET input file was successful" << endl;
	}
	else {
		cout << "Processing EPANET input file was not successful" << endl;
		return 1;
	}


	// Find dead end branches in the network
	N_DEbranches = DEFIND(&net);
	if (N_DEbranches) {
		cout << "Found " << N_DEbranches << " dead-end branches in the network --> See DEpipeID.txt" << endl;
	}
	else {
		cout << "Found no dead-end branches in the network" << endl;
		return 1;
	}


	// Find cross junctions in the network
	N_xjuncts = XJFIND(&net);
	cout << "Found " << N_xjuncts << " cross junctions in the network --> See XjuncIDs.txt" << endl;


	// Process EPANET report file
	error = OpenEPANETrpt(RPTfileName, &net);
	if (!error) {
		cout << "Processing EPANET report file was successful --> Starting Water Quality simulations" << endl;
	}
	else {
		cout << "Processing EPANET report file was not successful" << endl;
		return 1;
	}

	//  Run Water Quality simulations for dead-end branches
	error = WQSIM(WUDESIMinp, &net);
	if (!error) {
		cout << "Water Quality simulations finished successfuly" << endl;
	}
	else {
		cout << "Water Quality simulations were not successful" << endl;
		return 1;
	}


	return 0;
}