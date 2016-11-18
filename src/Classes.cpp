/*
Project:     WUDESIM ver. 1 BETA
File:        Classes.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: Class initialization with network defaults.
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

all_links::all_links() {
	node_1_conn = 0;	   //number of pipe connections to node 1
	node_2_conn = 0;	   //number of pipe connections to node 2
};


all_options::all_options() {

	Rel_Diffusivity = 1;
	Rel_Viscosity = 1;
	Flow_UNITS = "GPM";
	QUAL_TAG = "NONE";
	QUAL_UNIT = "";

	Flow_unit_conv = 0.0000630901964;
	unit_sys = 0;  //0-US 1-SI
};

all_times::all_times() {
	Duration_hr = 0;
	Duration_min = 0;
	Hyd_step_hr = 1;
	Hyd_step_min = 0;
	Qual_step_hr = 0;
	Qual_step_min = 0;
	Rep_step_hr = 1;
	Rep_step_min = 0;
	Rep_start_hr = 0;
	Rep_start_min = 0;
	N_steps;
};

all_reactions::all_reactions() {
	Bulk_coeff = 0.;
	Wall_coeff = 0.;
	Bulk_order = 1.;
	Wall_order = 1.;
	Lim_pot = 0.;
};

dead_end_branch::dead_end_branch() {
	branch_size = 0;						//Branch size;
};
