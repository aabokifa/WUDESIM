/*
Project:     WUDESIM ver. 1 BETA
File:        OpenEPANETrpt.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
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


#include "WUDESIM.h"
#include "Classes.h"
#include "WUDESIM_CORE.h"
#include "Utilities.h"
#include "WRITING_FUN.h"
#include "epanet2.h"


using namespace std;

int RUN_EPANET_SIM(Network* net) {

	const char* INPfileName = net->EPANET_INP;
	const char* RPTfileName = net->EPANET_RPT;

	// Initialize vectors for dead-end branches data
	int N_branches  = net->DE_branches.size();	
	
	for (int branch = 0;branch < N_branches;branch++) {

		int n_rows    = net->DE_branches[branch].branch_size;
		int n_columns = net->times.N_steps;

		net->DE_branches[branch].pipe_flow_EPANET.resize(n_rows, vector<double>(n_columns, 0.));
		net->DE_branches[branch].boundary_C_EPANET.resize(n_rows, vector<double>(n_columns, 0.));
		net->DE_branches[branch].terminal_C_EPANET.resize(n_rows, vector<double>(n_columns, 0.));
		net->DE_branches[branch].terminal_id.resize(n_rows);
		net->DE_branches[branch].bound_id.resize(n_rows);
	}	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Open EPANET 
	ENopen(INPfileName, RPTfileName, (char*)"");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Extract link flow results
	int N_links;
	long t, tstep;
	float link_flow;
	char  link_id[16];
	int Tstep;
	
	int second_interval = (net->times.Hyd_step_hr * 60 + net->times.Hyd_step_min) * 60;

	ENgetcount(EN_LINKCOUNT, &N_links);
	ENopenH();
	ENinitH(0);
	do {
		ENrunH(&t);


		if (t % second_interval == 0) {
			
			Tstep = t / second_interval;

			for (int i = 1;i <= N_links;i++) {

				ENgetlinkvalue(i, EN_FLOW, &link_flow);
				ENgetlinkid(i, link_id);


				for (int branch = 0;branch < N_branches;branch++) {

					for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

						if (link_id == net->DE_branches[branch].pipe_id[pipe]) {
							//Read flow data
							net->DE_branches[branch].pipe_flow_EPANET[pipe][Tstep] = link_flow;
						}
					}
				}
			}
		}
		
		ENnextH(&tstep);
	} while (tstep > 0);
	ENcloseH();
	ENsolveH(); // A complete hydraulic simulation is required before running a WQ simulation

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Based on the flow in each pipe, determine which node is boundary (inlet) and which is terminal (outlet)

	for (int branch = 0;branch < N_branches;branch++) {

		for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

			string bound_node, term_node;

			if (net->DE_branches[branch].pipe_flow_EPANET[pipe][Tstep] > 0) {

				bound_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_1;
				term_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_2;

				net->DE_branches[branch].bound_id[pipe]    = bound_node;
				net->DE_branches[branch].terminal_id[pipe] = term_node;
			}

			else {

				bound_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_2;
				term_node = net->pipes[net->DE_branches[branch].pipe_index[pipe]].node_1;

				net->DE_branches[branch].bound_id[pipe]    = bound_node;
				net->DE_branches[branch].terminal_id[pipe] = term_node;

			}
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Extract node quality results
	int N_nodes;
	long q, qstep;
	float node_qual;
	char node_id[16];

	ENgetcount(EN_NODECOUNT, &N_nodes);
	ENopenQ();
	ENinitQ(0);
	do {
		
		ENrunQ(&q);

		if (q % second_interval == 0) {

			Tstep = q / second_interval;

			for (int i = 1;i < N_nodes;i++) {
				ENgetnodevalue(i, EN_QUALITY, &node_qual);
				ENgetnodeid(i, node_id);


				for (int branch = 0;branch < N_branches;branch++) {

					for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

						if (node_id == net->DE_branches[branch].bound_id[pipe]) { net->DE_branches[branch].boundary_C_EPANET[pipe][Tstep] = node_qual; }
						if (node_id == net->DE_branches[branch].terminal_id[pipe]) { net->DE_branches[branch].terminal_C_EPANET[pipe][Tstep] = node_qual; }
					}

				}
			}


		}

		ENnextQ(&qstep);

	} while (qstep > 0);
	ENcloseQ();	
	ENsolveQ();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Save report file and close EPANET
	ENreport();
	ENclose();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Convert Flow Units

	for (int r_step = 0;r_step < net->times.N_steps;++r_step) {

		for (int branch = 0;branch < N_branches;branch++) {

			for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {
				//Convert flow units
				net->DE_branches[branch].pipe_flow_EPANET[pipe][r_step] *= net->options.Flow_unit_conv;

				//Correct negative flows
				if (net->DE_branches[branch].pipe_flow_EPANET[pipe][r_step] < 0) { net->DE_branches[branch].pipe_flow_EPANET[pipe][r_step] *= -1; }
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Save dead-end length and diameter 

	for (int branch = 0;branch < N_branches;branch++) {
		
		net->DE_branches[branch].length.resize(net->DE_branches[branch].branch_size);
		net->DE_branches[branch].diameter.resize(net->DE_branches[branch].branch_size);

		for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;++pipe) {

			net->DE_branches[branch].length[pipe]   = net->pipes[net->DE_branches[branch].pipe_index[pipe]].length;
			net->DE_branches[branch].diameter[pipe] = net->pipes[net->DE_branches[branch].pipe_index[pipe]].diameter;

		}
	}
	return 0;
}