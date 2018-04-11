/*
Project:     WUDESIM ver. 1 BETA
File:        XJFIND.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: This function studies the connectivity of each junction in the network to determine if it has more than four pipe connections,
             and generates a list of cross junction IDs to the file “XjuncIDs.txt”.
*/
#include <iostream> 
#include <fstream>  
#include <vector>  
#include <map>
#include <algorithm>
#include <sstream>  
#include <string>  
#include <iterator>
#include <numeric>


#include "Classes.h"
#include "WUDESIMmain.h"
#include "Utilities.h"

using namespace std;


int XJFIND(Network* net)
{
	//Find X-junctions in the network
	for (int i = 0;i < net->junctions.size();i++) {

		int XjuncCond = 0; //Condition for a X junction = must have more than four pipe connections

		for (int j = 0;j < net->pipes.size();j++) {

			if (net->junctions[i].id == net->pipes[j].node_1 || net->junctions[i].id == net->pipes[j].node_2) {
				XjuncCond++;
			}

		}

		if (XjuncCond >= 4) { net->XJunctions.push_back(net->junctions[i].id); }
	}
	

	// Write X-junction ids to XjuncIDs.txt
	ofstream ofs;
	ofs.open("XjuncIDs.txt", ios::out | ios::trunc);
	if (ofs.is_open()) {
		ofs << "Cross junction IDs:" << endl;
		for (int j = 0;j <  net->XJunctions.size();++j) {
			ofs << net->XJunctions[j] << endl;
		}
		ofs.close();
	}

	return net->XJunctions.size();

}