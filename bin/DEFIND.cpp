/*
File:   DEFIND.cpp
Author: Ahmed Abokifa
Date:   10/25/2016
Desc:	This is the entry point to the software. It starts by importing the EPANET input file specified by the user in WUDESIM.inp. 
		It then reads all the network elements data for pipes, junctions, tanks, reservoirs, pumps, and valves. It also reads all 
		the simulation parameters on time steps, reactions, and other simulation options. The code then studies the connectivity 
		of each junction to determine if it is a dead end, and generates a list of dead end pipe IDs is to the file “DEpipeID.txt”. 
		The code also identifies all cross-junctions in the network, and writes a list of their IDs to the file “XjuncIDs.txt”. 
		The code then imports the flow rates and boundary conditions for all dead ends in the network from the EPANET output report
		file specified by the user in "WUDESIM.inp". It then calls the water quality simulation module located in the source file 
		"WQSIM.cpp". The data is propagated to the function through temporary binary files that are deleted right after the simulation
		is finished.
*/


//#include "stdafx.h" 
#include <iostream> 
#include <fstream>  
#include <vector>  
#include <map>
#include <algorithm>
#include <sstream>  
#include <string>  
#include <iterator>
#include <numeric>
#include "Headers.h"
#include <stdio.h>
#include <tchar.h>
#include <SDKDDKVer.h>

using namespace std;

ofstream logfile;
ofstream myrptfile;
ofstream ofs;


class all_links {
public:
	string id;			    //Pipe ID
	double length;			//Pipe length
	double diameter;		//Pipe diameter
	string node_1;			//node 1 of pipe
	string node_2;			//node 2 of pipe
	int node_1_conn=0;		//number of pipe connections to node 1
	int node_2_conn=0;		//number of pipe connections to node 2
};

class dead_end_branch {
public:	
	int branch_size=0;							//Branch size;
	vector<string> pipe_id;						//IDs of the pipes in the DE branch
	vector<string> terminal_id;					//IDs of the terminal junctions in the DE branch
	vector<int> pipe_index;					    //Indices of the pipes in the DE branch
	vector<vector<double>> pipe_flow;			//Flow profile of dead end pipes (row=pipe/ column=flow@time)
	vector<vector<double>> boundary;			//Boundary condition profile     (row=pipe /column=concentration@time)
	vector<vector<double>> terminal;            //terminal concentration profile (row=pipe /column=concentration@time)
};

int main()
{

	// Write WUDESIM.log file
	logfile.open("WUDESIM.log", ios::out | ios::trunc);
	logfile << "---------WUDESIM log-----------" << endl;

	//Open WUDESIM.rpt file
	myrptfile.open("WUDESIM.rpt", ios::out | ios::trunc);
	myrptfile << "---------WUDESIM Results-----------" << endl;

	// Import WUDESIM.inp file
	vector<string> WUDESIMinp;
	WUDESIMinp = ImportFile("WUDESIM.inp");
	vector<string> Headers1 = { "[INP_FILE]","[RPT_FILE]","[CORR_FACTS]","[STOC_DEMANDS]","[END]" };

	vector<int> index1;
	for (int j = 0;j < Headers1.size();++j) {
		for (int i = 0;i < WUDESIMinp.size();++i) { if (find(Headers1[j], WUDESIMinp[i])) { index1.push_back(i); } };
	}
	sort(index1.begin(), index1.end());

	string INPfileName;
	string RPTfileName;
	vector<string> dum;   //dummy string variable
	dum = InputData(index1, WUDESIMinp, "[INP_FILE]");
	INPfileName = dum[0];
	dum = InputData(index1, WUDESIMinp, "[RPT_FILE]");
	RPTfileName = dum[0];

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Import EPANET input file (.inp)
	logfile << " Importing EPANET input file data" << endl << endl;
	vector<string> A1;
	A1 = ImportFile(INPfileName);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Find section locations

	vector<string> Headers = {
		"[TITLE]","[JUNCTIONS]","[RESERVOIRS]","[TANKS]","[PIPES]","[PUMPS]","[VALVES]","[EMITTERS]",
		"[CURVES]","[PATTERNS]","[ENERGY]","[STATUS]","[CONTROLS]","[RULES]","[DEMANDS]",
		"[QUALITY]","[REACTIONS]","[SOURCES]","[MIXING]",
		"[OPTIONS]","[TIMES]","[REPORT]",
		"[COORDINATES]","[VERTICES]","[LABELS]","[BACKDROP]","[TAGS]",
		"[END]" };

	vector<int> index;
	for (int j = 0;j < Headers.size();++j) {
		for (int i = 0;i < A1.size();++i) { if (find(Headers[j], A1[i])) { index.push_back(i); } };
	}
	sort(index.begin(), index.end());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Pipes data from input file

	vector<string> pipe_data;
	pipe_data = InputData(index, A1, "[PIPES]");

	int N_pipes = pipe_data.size();    //Number of pipes

	vector<all_links> pipes(N_pipes);

	for (int i = 0;i < N_pipes;++i) {
		istringstream iss(pipe_data[i]);
		iss >> pipes[i].id >> pipes[i].node_1 >> pipes[i].node_2 >> pipes[i].length >> pipes[i].diameter;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Junctions data from input file

	vector<string> node_data;
	node_data = InputData(index, A1, "[JUNCTIONS]");

	int N_nodes = node_data.size();

	vector<string> node_id(N_nodes);
	vector<double> node_elev(N_nodes);
	vector<double> node_demand(N_nodes);
	vector<string> source_nodes;
	int N_sources = 0;

	for (int i = 0;i < N_nodes;++i) {
		istringstream iss(node_data[i]);
		iss >> node_id[i] >> node_elev[i] >> node_demand[i];
		if (node_demand[i] < 0) { source_nodes.push_back(node_id[i]); ++N_sources; }
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Read Tanks data from input file

	vector<string> tanks_data;
	tanks_data = InputData(index, A1, "[TANKS]");

	int N_tanks = tanks_data.size();

	vector<string> tank_id(N_tanks);
	for (int i = 0;i < N_tanks;++i) {
		istringstream iss(tanks_data[i]);
		iss >> tank_id[i];
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Read Reservoirs data from input file
	vector<string> reserv_data;
	reserv_data = InputData(index, A1, "[RESERVOIRS]");

	int N_reserv = reserv_data.size();

	vector<string> reserv_id(N_reserv);
	for (int i = 0;i < N_reserv;++i) {
		istringstream iss(reserv_data[i]);
		iss >> reserv_id[i];
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Read Pumps data from input file

	vector<string> pumps_data;
	pumps_data = InputData(index, A1, "[PUMPS]");

	int N_pumps = pumps_data.size();

	vector<string> pump_id(N_pumps);
	vector<string> pump_start(N_pumps);
	vector<string> pump_end(N_pumps);

	for (int i = 0;i < N_pumps;++i) {
		istringstream iss(pumps_data[i]);
		iss >> pump_id[i] >> pump_start[i] >> pump_end[i];
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Read Valves data from input file

	vector<string> valves_data;
	valves_data = InputData(index, A1, "[VALVES]");

	int N_valves = valves_data.size();

	vector<string> valve_id(N_valves);
	vector<string> valve_start(N_valves);
	vector<string> valve_end(N_valves);

	for (int i = 0;i < N_valves;++i) {
		istringstream iss(valves_data[i]);
		iss >> valve_id[i] >> valve_start[i] >> valve_end[i];
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Find number of connections for node_1 and node_2 of each pipe
	vector<string> node_1_2(2 * N_pipes);
	int j = 0;
	for (int i = 0;i < N_pipes;++i) {
		node_1_2[j] = pipes[i].node_1;
		j++;
		node_1_2[j] = pipes[i].node_2;
		j++;
	}

	for (int i = 0;i < N_pipes;i++) {

		//Find number of connections for node_1 and node_2 of each pipe
		for (int j = 0;j < node_1_2.size();++j) {
			if (compare(pipes[i].node_1, node_1_2[j])) {
				pipes[i].node_1_conn++;
			}
			if (compare(pipes[i].node_2, node_1_2[j])) {
				pipes[i].node_2_conn++;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Find dead end branches

	vector<dead_end_branch> DE_branches;

	int branch = 0;        // First DE branch
	int N_branches = 0;

	string dead_node, other_node;

	int other_node_conn;

	for (int i = 0;i < N_pipes;i++) {

		// Condition 1: Check if one of the two nodes is connected to only one pipe
		if (pipes[i].node_1_conn == 1) {

			dead_node = pipes[i].node_1;
			other_node = pipes[i].node_2;
			other_node_conn = pipes[i].node_2_conn; //Number of connections to the other node
		}
		else if (pipes[i].node_2_conn == 1) {

			dead_node = pipes[i].node_2;
			other_node = pipes[i].node_1;
			other_node_conn = pipes[i].node_1_conn; //Number of connections to the other node
		}

		else { goto nextBranch; } //if it's not a dead-end then go directly to the next pipe


		//Condition 2: Check DE node is not a source
		for (int j = 0;j < N_sources;++j) {
			if (compare(dead_node, source_nodes[j])) { goto nextBranch; }
		}

		// Condition 3: Check DE node is not a tank
		for (int j = 0;j < N_tanks;++j) {
			if (compare(dead_node, tank_id[j])) { goto nextBranch; }
		}

		// Condition 4: Check DE node is not a reservoir
		for (int j = 0;j < N_reserv;++j) {
			if (compare(dead_node, reserv_id[j])) { goto nextBranch; }
		}

		// Condition 5: Check DE node is not a pump
		for (int j = 0;j < N_pumps;++j) {
			if (compare(dead_node, pump_start[j])) { goto nextBranch; }
			if (compare(dead_node, pump_end[j])) { goto nextBranch; }
		}

		// Condition 6: Check DE node is not a valve
		for (int j = 0;j < N_valves;++j) {
			if (compare(dead_node, valve_start[j])) { goto nextBranch; }
			if (compare(dead_node, valve_end[j])) { goto nextBranch; }
		}

		// if you come to here then you certainly are dead end!
		N_branches++;
		DE_branches.resize(N_branches);
		DE_branches[branch].pipe_id.push_back(pipes[i].id);
		DE_branches[branch].pipe_index.push_back(i);
		DE_branches[branch].branch_size++;

		//Now check to see if preceding pipes are to be added to the branch
		if (other_node_conn == 2) {

			int k = i;

		anotherPipe:; //Will come back to here if more pipes are to be added

			for (int p = 0;p < N_pipes;p++) {

				if (compare(pipes[p].node_1, other_node) && p != k) {

					DE_branches[branch].pipe_id.push_back(pipes[p].id);
					DE_branches[branch].pipe_index.push_back(p);
					DE_branches[branch].branch_size++;

					if (pipes[p].node_2_conn == 2) {
						other_node = pipes[p].node_2;
						k = p;
						goto anotherPipe;
					}
					else { branch++; goto nextBranch; }
				}

				else if (compare(pipes[p].node_2, other_node) && p != k) {

					DE_branches[branch].pipe_id.push_back(pipes[p].id);
					DE_branches[branch].pipe_index.push_back(p);
					DE_branches[branch].branch_size++;


					if (pipes[p].node_1_conn == 2) {
						other_node = pipes[p].node_1;
						k = p;
						goto anotherPipe;
					}
					else { branch++; goto nextBranch; }


				}
			}

		}
		else { branch++; goto nextBranch; }

	nextBranch:;
	}

	logfile << "Found " << N_branches << " Dead End Branches  in the network!" << endl << endl;


	// Write dead end branch ids
	ofs.open("DEpipeID.txt", ios::out | ios::trunc);
	if (ofs.is_open()) {
		ofs << "Deadend pipe IDs:" << endl;

		for (int branch = 0;branch < N_branches;branch++) {
			ofs << "DE Branch No. " << branch + 1 << " :" << endl;

			for (int pipe = 0;pipe < DE_branches[branch].branch_size;pipe++) {
				ofs << DE_branches[branch].pipe_id[pipe] << endl;

			}
		}
		ofs.close();
	}
	logfile << "Dead-End IDs written to 'DEpipeID.txt' " << endl << endl;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Find X-junctions in the network

	logfile << "Finding X-junctions in the network" << endl << endl;

	vector<string> XjuncID;
	vector<int> XjuncCond(node_id.size());

	for (int i = 0;i < node_id.size();++i) {
		for (int j = 0;j < node_1_2.size();++j) {
			if (compare(node_id[i], node_1_2[j])) {
				XjuncCond[i]++;
			}
		}
		if (XjuncCond[i] >= 4) { XjuncID.push_back(node_id[i]); }
	}
	logfile << "Found " << XjuncID.size() << " X-junctions in the network!" << endl << endl;

	// Write X-junction ids
	ofs.open("XjuncIDs.txt", ios::out | ios::trunc);
	if (ofs.is_open()) {
		ofs << "Cross junction IDs:" << endl;
		for (int j = 0;j < XjuncID.size();++j) {
			ofs << XjuncID[j] << endl;
		}
		ofs.close();
	}
	logfile << "X-junction IDs written to 'XjuncIDs.txt' " << endl << endl;


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Options data from input file
	string dummy;

	vector<string> options_data;
	options_data = InputData(index, A1, "[OPTIONS]");

	double Rel_Diffusivity;
	double Rel_Viscosity;
	string Flow_UNITS;
	string QUAL_TAG;

	for (int i = 0;i < options_data.size();++i) {
		istringstream iss(options_data[i]);

		if (find("Viscosity", options_data[i])) { iss >> dummy >> Rel_Viscosity; }
		if (find("Diffusivity", options_data[i])) { iss >> dummy >> Rel_Diffusivity; }
		if (find("UNITS", options_data[i])) { iss >> dummy >> Flow_UNITS; }
		if (find("QUALITY", options_data[i])) { iss >> dummy >> QUAL_TAG; }
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Flow Unit Conversion to m3/sec
	double Flow_unit_conv;
	int unit_sys; //0-US 1-SI

	if (Flow_UNITS == "CFS") { Flow_unit_conv = 0.028316847; unit_sys = 0; }
	else if (Flow_UNITS == "GPM") { Flow_unit_conv = 0.0000630901964; unit_sys = 0; }
	else if (Flow_UNITS == "MGD") { Flow_unit_conv = 0.043812636574; unit_sys = 0; }
	else if (Flow_UNITS == "IMGD") { Flow_unit_conv = 0.052616782407; unit_sys = 0; }
	else if (Flow_UNITS == "AFD") { Flow_unit_conv = 0.014276410185; unit_sys = 0; }
	else if (Flow_UNITS == "LPS") { Flow_unit_conv = 0.001; unit_sys = 1; }
	else if (Flow_UNITS == "LPM") { Flow_unit_conv = 0.000016666666667; unit_sys = 1; }
	else if (Flow_UNITS == "MLD") { Flow_unit_conv = 0.0115740741; unit_sys = 1; }
	else if (Flow_UNITS == "CMH") { Flow_unit_conv = 0.00027777777778; unit_sys = 1; }
	else if (Flow_UNITS == "CMD") { Flow_unit_conv = 0.000011574074074; unit_sys = 1; }

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Times data from input file

	vector<string> times_data;
	times_data = InputData(index, A1, "[TIMES]");

	char dummy1;
	int Duration_hr;
	int Duration_min;
	int Hyd_step_hr;
	int Hyd_step_min;
	int Qual_step_hr;
	int Qual_step_min;
	int Rep_step_hr;
	int Rep_step_min;
	int Rep_start_hr;
	int Rep_start_min;

	for (int i = 0;i < times_data.size();++i) {
		istringstream iss(times_data[i]);

		// Read duration
		if (find("Duration", times_data[i])) { iss >> dummy >> Duration_hr >> dummy1 >> Duration_min; goto nexttime; }

		// Read Hydraulic Time step
		if (find("Hydraulic Timestep", times_data[i])) { iss >> dummy >> dummy >> Hyd_step_hr >> dummy1 >> Hyd_step_min; goto nexttime; }

		// Read Quality Time step
		if (find("Quality Timestep", times_data[i])) { iss >> dummy >> dummy >> Qual_step_hr >> dummy1 >> Qual_step_min;goto nexttime; }

		// Read Report Time step
		if (find("Report Timestep", times_data[i])) { iss >> dummy >> dummy >> Rep_step_hr >> dummy1 >> Rep_step_min;goto nexttime; }

		// Read Report Start
		if (find("Report Start", times_data[i])) { iss >> dummy >> dummy >> Rep_start_hr >> dummy1 >> Rep_start_min;goto nexttime; }
	nexttime:;
	}

	int N_steps = ((Duration_hr - Rep_start_hr) * 60 + (Duration_min - Rep_start_min)) / (Rep_step_hr * 60 + Rep_step_min) + 1;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Read Global Reactions data from input file

	vector<string> reactions_data;
	reactions_data = InputData(index, A1, "[REACTIONS]");

	double Bulk_coeff = 0;
	double Wall_coeff = 0;
	double Bulk_order = 1;
	double Wall_order = 1;
	double Lim_pot = 0;

	for (int i = 0;i < reactions_data.size();++i) {
		istringstream iss(reactions_data[i]);

		// Read duration
		if (find("Global Bulk", reactions_data[i])) { iss >> dummy >> dummy >> Bulk_coeff; }  // Read Bulk Coeff (/day)
		if (find("Global Wall", reactions_data[i])) { iss >> dummy >> dummy >> Wall_coeff; }  // Read Wall Coeff (length/day)
		if (find("Order Bulk", reactions_data[i])) { iss >> dummy >> dummy >> Bulk_order; }
		if (find("Order Wall", reactions_data[i])) { iss >> dummy >> dummy >> Wall_order; }
		if (find("Limiting Potential", reactions_data[i])) { iss >> dummy >> dummy >> Lim_pot; }

	}

	if (find("NONE", QUAL_TAG)) { Bulk_coeff = 0; Wall_coeff = 0; }
	if (find("AGE", QUAL_TAG)) { Bulk_coeff = 1*24*3600; Bulk_order = 0; Wall_coeff = 0; }
	if (find("TRACE", QUAL_TAG)) { logfile << "WUDESIM does not handle Trace quality simulations" << endl; exit(EXIT_FAILURE); }

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Import EPANET report file (.rpt)

	logfile << " Importing EPANET report file data" << endl << endl;

	vector<string> A2;
	A2 = ImportFile(RPTfileName);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Find header locations

	vector<int> k_link;
	vector<int> k_node;

	int position = 0; //To continue search from last position
	for (int r_step = 0;r_step < N_steps;++r_step) {

		double Time = (r_step)*(Rep_step_hr * 60 + Rep_step_min) + (Rep_start_hr * 60 + Rep_start_min);
		int Time_hr = floor(Time / 60.);
		int Time_min = Time - Time_hr * 60;

		string Time_min_str;
		if (Time_min > 9) { Time_min_str = to_string(Time_min); }
		else { Time_min_str = "0" + to_string(Time_min); }

		string Time_hr_str = to_string(Time_hr);

		string link_str = "Link Results at " + Time_hr_str + ":" + Time_min_str + " Hrs:";
		string node_str = "Node Results at " + Time_hr_str + ":" + Time_min_str + " Hrs:";


		for (int i = position;i < A2.size();++i) {
			++position;
			if (find(node_str, A2[i])) { k_node.push_back(i); goto find_link; }
		}
	find_link:;
		for (int i = position;i < A2.size();++i) {
			++position;
			if (find(link_str, A2[i])) { k_link.push_back(i); goto find_node; }
		}
	find_node:;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Flow Data for dead end branches

	for (int branch = 0;branch < N_branches;branch++) {

		int n_rows = DE_branches[branch].branch_size;
		int n_columns = N_steps;

		DE_branches[branch].pipe_flow.resize(n_rows, vector<double>(n_columns, 0.));
		DE_branches[branch].boundary.resize(n_rows, vector<double>(n_columns, 0.));
		DE_branches[branch].terminal.resize(n_rows, vector<double>(n_columns, 0.));
		DE_branches[branch].terminal_id.resize(n_rows);
	}

	for (int r_step = 0;r_step < N_steps;++r_step) {

		int start_line = k_link[r_step];
		int end_line = (r_step < N_steps - 1) ? k_node[r_step + 1] : A2.size();

		for (int i = start_line;i < end_line;++i) {
			istringstream iss(A2[i]);
			string id;
			iss >> id;

			for (int branch = 0;branch < N_branches;branch++) {

				for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {

					if (id == DE_branches[branch].pipe_id[pipe]) {

						//Read flow data
						iss >> DE_branches[branch].pipe_flow[pipe][r_step];

					}
				}

			}
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Read Boundary & Terminal Concentration

	for (int r_step = 0;r_step < N_steps;++r_step) {

		int start_line = k_node[r_step];
		int end_line = k_link[r_step];

		for (int branch = 0;branch < N_branches;branch++) {

			for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {

				for (int i = start_line;i < end_line;++i) {

					istringstream iss(A2[i]);
					string id;
					iss >> id; //Read node id;

					string bound_node, term_node;

					if (DE_branches[branch].pipe_flow[pipe][r_step] > 0) {

						bound_node = pipes[DE_branches[branch].pipe_index[pipe]].node_1;
						term_node = pipes[DE_branches[branch].pipe_index[pipe]].node_2;
						DE_branches[branch].terminal_id[pipe] = term_node;
					}

					else {

						bound_node = pipes[DE_branches[branch].pipe_index[pipe]].node_2;
						term_node = pipes[DE_branches[branch].pipe_index[pipe]].node_1;
						DE_branches[branch].terminal_id[pipe] = term_node;

					}

					double dum;
					if (id == bound_node) { iss >> dum >> dum >> dum >> DE_branches[branch].boundary[pipe][r_step]; }
					if (id == term_node) { iss >> dum >> dum >> dum >> DE_branches[branch].terminal[pipe][r_step]; }
				}
			}
		}

	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int r_step = 0;r_step < N_steps;++r_step) {

		for (int branch = 0;branch < N_branches;branch++) {

			for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {
				//Convert flow units
				DE_branches[branch].pipe_flow[pipe][r_step] *= Flow_unit_conv;

				//Correct negative flows
				if (DE_branches[branch].pipe_flow[pipe][r_step] < 0) { DE_branches[branch].pipe_flow[pipe][r_step] *= -1; }
			}
		}
	}

	// Pipe length and diameter conversion
	for (int i = 0;i < N_pipes;++i) {
		if (unit_sys == 0) {
			pipes[i].length *= 0.3048;     //ft-->m
			pipes[i].diameter *= 0.0254;   //in-->m
		}
		else {
			pipes[i].length *= 1;         //m-->m
			pipes[i].diameter *= 0.001;   //mm-->m
		}
	}

	// Reaction Unit Conversion
	// Bulk
	Bulk_coeff *= 1 / (24 * 3600); // 1/day --> 1/sec
	// Wall
	if (unit_sys == 0) {
		if (Wall_order == 1) { Wall_coeff *= 0.3048 / (24 * 3600); }                       // ft/day --> m/sec
		else if (Wall_order == 0) { Wall_coeff /= (pow(0.3048, 2) * 24 * 3600); }        // 1/ft2/day --> 1/m2/sec
	}
	else if (unit_sys == 1) { Wall_coeff *= 1 / (24 * 3600); }  // m/day --> m/sec

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Write simulation info data to binary file
	logfile << "Writing input data for water quality simulation" << endl << endl;

	vector<double> SIMinfo(10);
	SIMinfo[0] = N_steps;								 //Number of steps (hrs);
	SIMinfo[1] = Rep_step_hr + Rep_step_min / 60.;		 //Report time step (hr)
	SIMinfo[2] = Qual_step_hr*3600. + Qual_step_min*60.; //Quality time step (sec)
	SIMinfo[3] = Bulk_coeff;                             //bulk decay coeff (/sec)
	SIMinfo[4] = Wall_coeff;                             //wall decay coefficient (m/sec) or (1/m2/sec)
	SIMinfo[5] = Rel_Diffusivity*1.2E-9;                 //Actual Diffusivity (m2/sec)
	SIMinfo[6] = Rel_Viscosity*1E-6;                     //Actual Kinematic Viscosity (m2/sec)
	SIMinfo[7] = Bulk_order;							 //Order of bulk reaction
	SIMinfo[8] = Wall_order;							 //Order of wall reaction
	SIMinfo[9] = Lim_pot;								 //Limiting concentration potential

	ofs.open("SIMinfo.bin", ios::out | ios::binary | ios::trunc);
	if (ofs.is_open()) {
		ofs.write(reinterpret_cast<char*>(SIMinfo.data()), (SIMinfo.size()) * sizeof(double));
		ofs.close();
	}
	else {
		logfile << "Error opening SIMinfo.bin file" << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	logfile << "Start water quality simulations" << endl << endl;


	// Write dead end information to binary files

	for (int branch = 0;branch < N_branches;branch++) {

		// Write Dead End branch size, pipe lengths, and diameters
		ofs.open("DEinfo.bin", ios::out | ios::binary | ios::trunc);

		if (ofs.is_open()) {

			ofs.write(reinterpret_cast<char*>(&DE_branches[branch].branch_size), sizeof(int));     //First write number of pipes in DE branch

			for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {

				ofs.write(reinterpret_cast<char*>(&pipes[DE_branches[branch].pipe_index[pipe]].length), sizeof(double));
				ofs.write(reinterpret_cast<char*>(&pipes[DE_branches[branch].pipe_index[pipe]].diameter), sizeof(double));

			}

			ofs.close();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// Write dead end flow data to binary file

		ofs.open("DEflow.bin", ios::out | ios::binary | ios::trunc);
		if (ofs.is_open()) {
			for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {
				for (int r_step = 0;r_step < N_steps;++r_step) {
					ofs.write(reinterpret_cast<char*>(&DE_branches[branch].pipe_flow[pipe][r_step]), sizeof(double));
				}
			}
			ofs.close();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Write dead end boundary concentration data to binary file
		// We will only write boundary concentration for the first pipe in the dead end branch

		ofs.open("DEboundary.bin", ios::out | ios::binary | ios::trunc);
		if (ofs.is_open()) {
			int first_pipe = DE_branches[branch].branch_size - 1;
			for (int r_step = 0;r_step < N_steps;++r_step) {
				ofs.write(reinterpret_cast<char*>(&DE_branches[branch].boundary[first_pipe][r_step]), sizeof(double));
			}

			ofs.close();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Write dead end terminal concentrations data to binary file

		ofs.open("DEterminal.bin", ios::out | ios::binary | ios::trunc);
		if (ofs.is_open()) {
			for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {
				for (int r_step = 0;r_step < N_steps;++r_step) {
					ofs.write(reinterpret_cast<char*>(&DE_branches[branch].terminal[pipe][r_step]), sizeof(double));
				}
			}
			ofs.close();
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Call Water Quality Simulation solver

		vector<string> term_id(DE_branches[branch].terminal_id.size());
		for (int pipe = 0;pipe < DE_branches[branch].branch_size;++pipe) {
			term_id[pipe]=DE_branches[branch].terminal_id[pipe];
		}

		WUDESIM(branch,term_id);

	}

	logfile << "Water Quality Simulations finished" << endl << endl;

	remove("DEboundary.bin");
	remove("DEterminal.bin");
	remove("DEflow.bin");
	remove("DEinfo.bin");
	remove("SIMinfo.bin");

	return 0;
}