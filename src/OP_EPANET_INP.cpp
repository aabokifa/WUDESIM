/*
Project:     WUDESIM ver. 1 BETA
File:        OpenEPANETinp.cpp
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

int OP_EPANET_INP(Network* net)
{

	// Import EPANET input file (.inp)
	string INPfileName = net->EPANET_INP;

	vector<string> EPANETinp;
	EPANETinp = ImportFile(INPfileName);
	
	if (EPANETinp.empty()) {
		WRITE_LOG_MSG(INPfileName + "is empty/corrupt!");
		return 1;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Find section locations

	vector<int> index;
	vector<string> Headers = {
		"[TITLE]","[JUNCTIONS]","[RESERVOIRS]","[TANKS]","[PIPES]","[PUMPS]","[VALVES]","[EMITTERS]",
		"[CURVES]","[PATTERNS]","[ENERGY]","[STATUS]","[CONTROLS]","[RULES]","[DEMANDS]",
		"[QUALITY]","[REACTIONS]","[SOURCES]","[MIXING]",
		"[OPTIONS]","[TIMES]","[REPORT]",
		"[COORDINATES]","[VERTICES]","[LABELS]","[BACKDROP]","[TAGS]",
		"[END]" };

	for (int j = 0;j < Headers.size();++j) {
		for (int i = 0;i < EPANETinp.size();++i) { if (find_str(Headers[j], EPANETinp[i])) { index.push_back(i); } };
	}
	sort(index.begin(), index.end());

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Pipes data from input file

	vector<string> pipe_data;
	pipe_data = InputData(index, EPANETinp, "[PIPES]");

	int N_pipes = pipe_data.size();    //Number of pipes
	
	if (N_pipes == 0) {
		return 1;
	}
	else {
		net->pipes.resize(N_pipes);

		for (int i = 0;i < N_pipes;++i) {
			istringstream iss(pipe_data[i]);
			iss >> net->pipes[i].id >> net->pipes[i].node_1 >> net->pipes[i].node_2 >> net->pipes[i].length >> net->pipes[i].diameter;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Junctions data from input file

	vector<string> node_data;
	node_data = InputData(index, EPANETinp, "[JUNCTIONS]");

	int N_nodes = node_data.size();
	
	if (N_nodes == 0) {
		return 1;
	}
	else {
		net->junctions.resize(N_nodes);

		for (int i = 0;i < N_nodes;++i) {
			istringstream iss(node_data[i]);
			iss >> net->junctions[i].id >> net->junctions[i].elev >> net->junctions[i].demand;			
		}
	}

	// Read supplemnt junction demands
	node_data = InputData(index, EPANETinp, "[DEMANDS]");
	if (node_data.size() != 0) {
		
		string junc_id;
		double junc_demand;
		for (int i = 0;i < node_data.size(); ++i) {
			istringstream iss(node_data[i]);
			iss >> junc_id >> junc_demand;
			for (int j = 0;j < N_nodes;++j) {
				if (compare_str(junc_id,net->junctions[j].id)) {
					net->junctions[i].demand = net->junctions[i].demand + junc_demand;
				}
			}
		}
	}
	
	for (int i = 0;i < N_nodes;++i) {

		if (net->junctions[i].demand < 0)
		{
			net->demand_sources.push_back(net->junctions[i].id);
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Read Tanks data from input file

	vector<string> tanks_data;
	tanks_data = InputData(index, EPANETinp, "[TANKS]");

	int N_tanks = tanks_data.size();

	net->tanks.resize(N_tanks);
	for (int i = 0;i < N_tanks;++i) {
		istringstream iss(tanks_data[i]);
		iss >> net->tanks[i];
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Read Reservoirs data from input file

	vector<string> reserv_data;
	reserv_data = InputData(index, EPANETinp, "[RESERVOIRS]");

	int N_reserv = reserv_data.size();

	net->reservoirs.resize(N_reserv);

	for (int i = 0;i < N_reserv;++i) {
		istringstream iss(reserv_data[i]);
		iss >> net->reservoirs[i];
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Read Pumps data from input file

	vector<string> pumps_data;
	pumps_data = InputData(index, EPANETinp, "[PUMPS]");

	int N_pumps = pumps_data.size();

	net->pumps.resize(N_pumps);

	for (int i = 0;i < N_pumps;++i) {
		istringstream iss(pumps_data[i]);
		iss >> net->pumps[i].id >> net->pumps[i].start >> net->pumps[i].end;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Read Valves data from input file

	vector<string> valves_data;
	valves_data = InputData(index, EPANETinp, "[VALVES]");

	int N_valves = valves_data.size();

	net->valves.resize(N_valves);

	for (int i = 0;i < N_valves;++i) {
		istringstream iss(valves_data[i]);
		iss >> net->valves[i].id >> net->valves[i].start >> net->valves[i].end;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Sources data from input file

	vector<string> sources_data;
	sources_data = InputData(index, EPANETinp, "[SOURCES]");

	int N_sources = sources_data.size();

	net->quality_sources.resize(N_sources);

	for (int i = 0;i < N_sources;++i) {
		istringstream iss(sources_data[i]);
		iss >> net->quality_sources[i];
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Options data from input file
	string dummy;

	vector<string> options_data;
	options_data = InputData(index, EPANETinp, "[OPTIONS]");

	for (int i = 0;i < options_data.size();++i) {
		istringstream iss(options_data[i]);

		if (find_str("Viscosity", options_data[i])) { iss >> dummy >> net->options.Rel_Viscosity; }
		if (find_str("Diffusivity", options_data[i])) { iss >> dummy >> net->options.Rel_Diffusivity; }
		if (find_str("UNITS", options_data[i])) { iss >> dummy >> net->options.Flow_UNITS; }
		if (find_str("QUALITY", options_data[i])) { iss >> dummy >> net->options.QUAL_TAG>> net->options.QUAL_UNIT; }
	}

	if (find_str("NONE", net->options.QUAL_TAG) || find_str("AGE", net->options.QUAL_TAG) || find_str("TRACE", net->options.QUAL_TAG)) {
		WRITE_LOG_MSG("WUDESIM can only take CHEMICAL water quality analysis");
		return 1;
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Times data from input file

	vector<string> times_data;
	times_data = InputData(index, EPANETinp, "[TIMES]");

	char dummy1;

	for (int i = 0;i < times_data.size();++i) {
		istringstream iss(times_data[i]);

		// Read duration
		if (find_str("Duration", times_data[i])) { iss >> dummy >> net->times.Duration_hr >> dummy1 >> net->times.Duration_min; goto nexttime; }
		
		// Read Hydraulic Time step
		if (find_str("Hydraulic Timestep", times_data[i])) { iss >> dummy >> dummy >> net->times.Hyd_step_hr >> dummy1 >> net->times.Hyd_step_min; goto nexttime; }

		// Read Quality Time step
		if (find_str("Quality Timestep", times_data[i])) { iss >> dummy >> dummy >> net->times.Qual_step_hr >> dummy1 >> net->times.Qual_step_min;goto nexttime; }
				
		// Read Report Time step
		if (find_str("Report Timestep", times_data[i])) { iss >> dummy >> dummy >> net->times.Rep_step_hr >> dummy1 >> net->times.Rep_step_min;goto nexttime; }

		// Read Report Start
		if (find_str("Report Start", times_data[i])) { iss >> dummy >> dummy >> net->times.Rep_start_hr >> dummy1 >> net->times.Rep_start_min;goto nexttime; }
	nexttime:;
	}

	// Check that it's not a single period snapshot simulation
	if (net->times.Duration_hr * 60 + net->times.Duration_min == 0) { WRITE_LOG_MSG("WUDESIM can't run single period snapshot analysis"); return 1; }
	
	// Check that the report time step and the hydraulic time step are equivalent
	if (net->times.Hyd_step_hr * 60 + net->times.Hyd_step_min != net->times.Rep_step_hr * 60 + net->times.Rep_step_min) {
		WRITE_LOG_MSG("The report time step must be equivalent to the hydraulic time step"); return 1;
	}

	// If no qualtiy step is defined use the default of one tenth of the hydraulic step
	if (net->times.Qual_step_hr + net->times.Qual_step_min == 0) {

		net->times.Qual_step_min = floor((net->times.Hyd_step_hr * 60 + net->times.Hyd_step_min) / 10);
	}
	
	// Get number of simulation steps (Note: both EPANET and WUDESIM simulate the entire simulation duration)
	net->times.N_steps = (net->times.Duration_hr * 60 + net->times.Duration_min) / (net->times.Hyd_step_hr * 60 + net->times.Hyd_step_min) + 1;
	
	// Get number of report steps
	net->times.N_steps_rep = ((net->times.Duration_hr - net->times.Rep_start_hr) * 60 + (net->times.Duration_min - net->times.Rep_start_min)) / (net->times.Rep_step_hr * 60 + net->times.Rep_step_min) + 1;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Read Global Reactions data from input file
	vector<string> reactions_data;
	reactions_data = InputData(index, EPANETinp, "[REACTIONS]");

	for (int i = 0;i < reactions_data.size();++i) {
		istringstream iss(reactions_data[i]);

		if (find_str("Global Bulk", reactions_data[i]))        { iss >> dummy >> dummy >> net->reactions.Bulk_coeff; }  // Read Bulk Coeff (/day)
		if (find_str("Global Wall", reactions_data[i]))        { iss >> dummy >> dummy >> net->reactions.Wall_coeff; }  // Read Wall Coeff (length/day)
		if (find_str("Order Bulk", reactions_data[i]))         { iss >> dummy >> dummy >> net->reactions.Bulk_order; }
		if (find_str("Order Wall", reactions_data[i]))         { iss >> dummy >> dummy >> net->reactions.Wall_order; }
		if (find_str("Limiting Potential", reactions_data[i])) { iss >> dummy >> dummy >> net->reactions.Lim_pot; }
	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Flow Unit Conversion to m3/sec
	if      (net->options.Flow_UNITS == "CFS") { net->options.Flow_unit_conv = 0.028316847;          net->options.unit_sys = 0; }
	else if (net->options.Flow_UNITS == "GPM") { net->options.Flow_unit_conv = 0.0000630901964;      net->options.unit_sys = 0; }
	else if (net->options.Flow_UNITS == "MGD") { net->options.Flow_unit_conv = 0.043812636574;       net->options.unit_sys = 0; }
	else if (net->options.Flow_UNITS == "IMGD"){ net->options.Flow_unit_conv = 0.052616782407;       net->options.unit_sys = 0; }
	else if (net->options.Flow_UNITS == "AFD") { net->options.Flow_unit_conv = 0.014276410185;       net->options.unit_sys = 0; }
	else if (net->options.Flow_UNITS == "LPS") { net->options.Flow_unit_conv = 0.001;                net->options.unit_sys = 1; }
	else if (net->options.Flow_UNITS == "LPM") { net->options.Flow_unit_conv = 0.000016666666667;    net->options.unit_sys = 1; }
	else if (net->options.Flow_UNITS == "MLD") { net->options.Flow_unit_conv = 0.0115740741;         net->options.unit_sys = 1; }
	else if (net->options.Flow_UNITS == "CMH") { net->options.Flow_unit_conv = 0.00027777777778;     net->options.unit_sys = 1; }
	else if (net->options.Flow_UNITS == "CMD") { net->options.Flow_unit_conv = 0.000011574074074;    net->options.unit_sys = 1; }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Pipe length and diameter conversion to meters
	for (int i = 0;i < net->pipes.size();++i) {
		if (net->options.unit_sys == 0) {
			net->pipes[i].length   *= 0.3048;     //ft-->m
			net->pipes[i].diameter *= 0.0254;     //in-->m
		}
		else {
			net->pipes[i].length   *= 1;          //m -->m
			net->pipes[i].diameter *= 0.001;      //mm-->m
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Reaction Unit Conversion
	// Bulk
	net->reactions.Bulk_coeff *= 1. / (24. * 3600.); // 1/day --> 1/sec
	
    // Wall
	if (net->options.unit_sys == 0) {
		if (net->reactions.Wall_order == 1) { net->reactions.Wall_coeff *= 0.3048 / (24. * 3600.); }                       // ft/day --> m/sec
		else if (net->reactions.Wall_order == 0) { net->reactions.Wall_coeff /= (pow(0.3048, 2) * 24. * 3600.); }          // 1/ft2/day --> 1/m2/sec
	}
	else if (net->options.unit_sys == 1) { net->reactions.Wall_coeff *= 1 / (24. * 3600.);}  // m/day --> m/sec

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;

}