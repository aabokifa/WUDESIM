/*
Project:     WUDESIM ver. 1 BETA
File:        OpenEPANETrpt.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: This function imports EPANET report file specified by the user, and reads flow rate, boundary condition, and terminal concentration
             profiles of all dead-end branches in the network.
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
#include <stdio.h>
#include <tchar.h>
#include <SDKDDKVer.h>

#include "Classes.h"
#include "WUDESIMmain.h"
#include "Utilities.h"

using namespace std;


int OpenEPANETrpt(string RPTfileName, Network* net) {	

	// Import EPANET report file (.rpt)

	vector<string> EPANETrpt;
	EPANETrpt = ImportFile(RPTfileName);

	if (EPANETrpt.empty()) {
		cout << RPTfileName << "is empty/corrupt!" << endl;
		return 1;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Find header locations

	vector<int> k_link;
	vector<int> k_node;

	int position = 0; //To continue search from last position
	for (int r_step = 0;r_step < net->times.N_steps;++r_step) {

		double Time = (r_step)*(net->times.Rep_step_hr * 60 + net->times.Rep_step_min) + (net->times.Rep_start_hr * 60 + net->times.Rep_start_min);
		int Time_hr = floor(Time / 60.);
		int Time_min = Time - Time_hr * 60;

		string Time_min_str;
		if (Time_min > 9) { Time_min_str = to_string(Time_min); }
		else { Time_min_str = "0" + to_string(Time_min); }

		string Time_hr_str = to_string(Time_hr);

		string link_str = "Link Results at " + Time_hr_str + ":" + Time_min_str + " Hrs:";
		string node_str = "Node Results at " + Time_hr_str + ":" + Time_min_str + " Hrs:";


		for (int i = position;i < EPANETrpt.size();++i) {
			++position;
			if (find(node_str, EPANETrpt[i])) { k_node.push_back(i); goto find_link; }
		}
	find_link:;
		for (int i = position;i < EPANETrpt.size();++i) {
			++position;
			if (find(link_str, EPANETrpt[i])) { k_link.push_back(i); goto find_node; }
		}
	find_node:;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Flow Data for dead end branches

	int N_branches = net->DE_branches.size();

	for (int branch = 0;branch < N_branches;branch++) {

		int n_rows = net->DE_branches[branch].branch_size;
		int n_columns = net->times.N_steps;

		net->DE_branches[branch].pipe_flow.resize(n_rows, vector<double>(n_columns, 0.));
		net->DE_branches[branch].boundary.resize(n_rows, vector<double>(n_columns, 0.));
		net->DE_branches[branch].terminal.resize(n_rows, vector<double>(n_columns, 0.));
		net->DE_branches[branch].terminal_id.resize(n_rows);
	}

	for (int r_step = 0;r_step < net->times.N_steps;++r_step) {

		int start_line = k_link[r_step];
		int end_line = (r_step < net->times.N_steps - 1) ? k_node[r_step + 1] : EPANETrpt.size();

		for (int i = start_line;i < end_line;++i) {
			istringstream iss(EPANETrpt[i]);
			string id;
			iss >> id;

			for (int branch = 0;branch < N_branches;branch++) {

				for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

					if (id == net->DE_branches[branch].pipe_id[pipe]) {

						//Read flow data
						iss >> net->DE_branches[branch].pipe_flow[pipe][r_step];
					}
				}
			}
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Boundary & Terminal Concentration

	for (int r_step = 0;r_step < net->times.N_steps;++r_step) {

		int start_line = k_node[r_step];
		int end_line = k_link[r_step];

		for (int branch = 0;branch < N_branches;branch++) {

			for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

				for (int i = start_line;i < end_line;++i) {

					istringstream iss(EPANETrpt[i]);
					string id;
					iss >> id; //Read node id;

					string bound_node, term_node;

					if (net->DE_branches[branch].pipe_flow[pipe][r_step] > 0) {

						bound_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_1;
						term_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_2;
						net->DE_branches[branch].terminal_id[pipe] = term_node;
					}

					else {

						bound_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_2;
						term_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_1;
						net->DE_branches[branch].terminal_id[pipe] = term_node;

					}

					double dum;
					if (id == bound_node) { iss >> dum >> dum >> dum >> net->DE_branches[branch].boundary[pipe][r_step]; }
					if (id == term_node) { iss >> dum >> dum >> dum >> net->DE_branches[branch].terminal[pipe][r_step]; }
				}
			}
		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Convert Flow Units

	for (int r_step = 0;r_step < net->times.N_steps;++r_step) {

		for (int branch = 0;branch < N_branches;branch++) {

			for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {
				//Convert flow units
				net->DE_branches[branch].pipe_flow[pipe][r_step] *= net->options.Flow_unit_conv;

				//Correct negative flows
				if (net->DE_branches[branch].pipe_flow[pipe][r_step] < 0) { net->DE_branches[branch].pipe_flow[pipe][r_step] *= -1; }
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Save dead-end length and diameter 

	for (int branch = 0;branch < N_branches;branch++) {
		net->DE_branches[branch].length.resize(net->DE_branches[branch].branch_size);
		net->DE_branches[branch].diameter.resize(net->DE_branches[branch].branch_size);

		for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

			net->DE_branches[branch].length[pipe]=	net->pipes[net->DE_branches[branch].pipe_index[pipe]].length;
			net->DE_branches[branch].diameter[pipe] = net->pipes[net->DE_branches[branch].pipe_index[pipe]].diameter;

		}
	}
	return 0;
}