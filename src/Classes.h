/*
Project:     WUDESIM ver. 1 BETA
File:        Classes.h
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: Header file for class definitions.
*/

#pragma once

#include <vector>  

using namespace std;

class all_links {
public:
	string id;			       //Pipe ID
	double length;			   //Pipe length
	double diameter;		   //Pipe diameter
	string node_1;			   //node 1 of pipe
	string node_2;			   //node 2 of pipe
	int node_1_conn = 0;	   //number of pipe connections to node 1
	int node_2_conn = 0;	   //number of pipe connections to node 2
};

class all_nodes {
public:
	string id;
	double elev;
	double demand;
};

class all_pumps {
public:
	string id;
	string start;
	string end;
};

class all_valves {
public:
	string id;
	string start;
	string end;
};

class all_options {
public:
	double Rel_Diffusivity=1;
	double Rel_Viscosity=1;
	string Flow_UNITS="GPM";
	string QUAL_TAG="NONE";
	string QUAL_UNIT="";

	double Flow_unit_conv = 0.0000630901964;
	int unit_sys = 0;  //0-US 1-SI
};

class all_times {
public:
	int Duration_hr=0;
	int Duration_min=0;
	int Hyd_step_hr=1;
	int Hyd_step_min=0;
	int Qual_step_hr=0;
	int Qual_step_min=0;
	int Rep_step_hr=1;
	int Rep_step_min=0;
	int Rep_start_hr=0;
	int Rep_start_min=0;
	int N_steps;
};

class all_reactions {
public:
	double Bulk_coeff = 0.;
	double Wall_coeff = 0.;
	double Bulk_order = 1;
	double Wall_order = 1;
	double Lim_pot = 0;
};

class dead_end_branch {
public:
	int branch_size = 0;						//Branch size;
	vector<string> pipe_id;						//IDs of the pipes in the DE branch
	vector<double> length;
	vector<double> diameter;
	vector<string> terminal_id;					//IDs of the terminal junctions in the DE branch
	vector<int> pipe_index;					    //Indices of the pipes in the DE branch
	vector<vector<double>> pipe_flow;			//Flow profile of dead end pipes (row=pipe/ column=flow@time)
	vector<vector<double>> boundary;			//Boundary condition profile     (row=pipe /column=concentration@time)
	vector<vector<double>> terminal;            //terminal concentration profile as simulated by EPANET
	vector<vector<double>> terminal_new;		//new terminal concentration as simulated by WUDESIM

};


class Network {
public:
	vector<all_links> pipes;

	vector<all_nodes> junctions;

	vector<string>demand_sources;

	vector<string>quality_sources;

	vector<string> tanks;

	vector<string> reservoirs;

	vector<all_pumps> pumps;

	vector<all_valves> valves;

	all_options options;

	all_times times;

	all_reactions reactions;

	vector<dead_end_branch> DE_branches;

	vector<string> XJunctions;

};
