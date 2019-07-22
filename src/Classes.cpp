/*
Project:     WUDESIM ver. 1 BETA
File:        Classes.cpp
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

using namespace std;

all_links::all_links() {
	node_1_conn = 0;	   //number of pipe connections to node 1
	node_2_conn = 0;	   //number of pipe connections to node 2
	diameter    = 0;
	length      = 0;
};

all_options::all_options() {

	Rel_Diffusivity = 1;
	Rel_Viscosity   = 1;
	Flow_UNITS      = "GPM";
	QUAL_TAG        = "NONE";
	QUAL_UNIT       = "";

	Flow_unit_conv = 0.0000630901964;
	unit_sys       = 0;  //0-US 1-SI
};

all_times::all_times() {
	Duration_hr   = 0;
	Duration_min  = 0;
	Hyd_step_hr   = 1;
	Hyd_step_min  = 0;
	Qual_step_hr  = 0;
	Qual_step_min = 0;
	Rep_step_hr   = 1;
	Rep_step_min  = 0;
	Rep_start_hr  = 0;
	Rep_start_min = 0;
	N_steps       = 0;
	N_steps_rep   = 0;
};

all_reactions::all_reactions() {
	Bulk_coeff = 0.;
	Wall_coeff = 0.;
	Bulk_order = 1.;
	Wall_order = 1.;
	Lim_pot    = 0.;
};

dead_end_branch::dead_end_branch() {
	branch_size = 0;						//Branch size;
	
};

dead_end_options::dead_end_options() {
	Corr_Fact_fl = 0; 
	
	conn_demand  = 0;
	seg_length   = 0;
	
	flow_corr_fl = 0;
	disp_corr_fl = 0;
	Rw_corr_fl = 0;

	Stoc_dem_fl     = 0;
	u1 = 0;
	u2 = 0;
	s1 = 0;
	s2 = 0;
	avg_int = 0;

	Dispersion_fl = 1;

}