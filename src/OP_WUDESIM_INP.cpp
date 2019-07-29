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

int OP_WUDESIM_INP(Network* net) {
	
	// Import WUDESIM.inp file
	string WUDESIMinpfileName = net->WUDESIM_INP;
	
	vector<string> WUDESIMinp;
	WUDESIMinp = ImportFile(WUDESIMinpfileName);

	if (WUDESIMinp.empty()) {
		WRITE_LOG_MSG(WUDESIMinpfileName + "is empty/corrupt!"); return 1;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Find section locations
	vector<int>    index;
	vector<string> Headers1 = { "[SOL_DISPERSION]","[CORR_FACTS]","[STOC_DEMANDS]","[WUDESIM_OUTPUT]","[END]" };

	for (int j = 0; j < Headers1.size(); ++j) {
		for (int i = 0; i < WUDESIMinp.size(); ++i) { if (find_str(Headers1[j], WUDESIMinp[i])) { index.push_back(i); } };
	}
	sort(index.begin(), index.end());

	if (index.size() < 5) { WRITE_LOG_MSG(WUDESIMinpfileName + " must have [SOL_DISPERSION], [CORR_FACTS], [STOC_DEMANDS], [WUDESIM_OUTPUT], and [END] sections"); return 1; }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Read dispersion data

	vector<string> dispersion_data;

	dispersion_data = InputData(index, WUDESIMinp, "[SOL_DISPERSION]");

	for (int i = 0; i < dispersion_data.size(); i++) {

		istringstream iss(dispersion_data[i]);
		string dummy;
		char Dispersion_fl;

		if (find_str("DISPERSION COEFF", dispersion_data[i])) {

			iss >> dummy >> dummy >> Dispersion_fl;			
			

			if (Dispersion_fl == 'N') {
				WRITE_LOG_MSG("	o	Solute  dispersion turned OFF ");
				net->DE_options.Dispersion_fl = 0;
			}
			else if (Dispersion_fl == 'T') {
				WRITE_LOG_MSG("	o	Solute  dispersion turned ON using Taylor's coefficients ");
				net->DE_options.Dispersion_fl = 1;
			}
			else if (Dispersion_fl == 'A') {
				WRITE_LOG_MSG("	o	Solute  dispersion turned ON using Lee 2004 average dispersion coefficients ");
				net->DE_options.Dispersion_fl = 2;
			}
			
			else {
				WRITE_LOG_MSG("	o	DISPERSION COEFF flag not recognized ");  return 1;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read correction factor data

	vector<string> Corr_fact_data;
	Corr_fact_data = InputData(index, WUDESIMinp, "[CORR_FACTS]");

	char Corr_fact_fl = 'N';
	
	double conn_demand;
	double seg_length;
	
	char flow_corr_fl = 'N';
	char disp_corr_fl = 'N';
	char Rw_corr_fl   = 'N';

	for (int i = 0; i < Corr_fact_data.size(); i++) {

		istringstream iss(Corr_fact_data[i]);
		string dummy;

		if (find_str("CORRECTION FACTORS", Corr_fact_data[i])) {

			iss >> dummy >> dummy >> Corr_fact_fl;

			if (Corr_fact_fl=='F') {       //Flow-based correction
				net->DE_options.Corr_Fact_fl = 1;
			}
			else if (Corr_fact_fl == 'S') { //Spacing-based correction
				net->DE_options.Corr_Fact_fl = 2;
			}
			else if (Corr_fact_fl == 'N') {  //No correction
				net->DE_options.Corr_Fact_fl = 0;
				WRITE_LOG_MSG("	o	Correction factors turned OFF ");
			}
			else {
				WRITE_LOG_MSG("	o	CORRECTION FACTORS flag not recognized ");  return 1;
			}
		}
		if (find_str("CONN DEM", Corr_fact_data[i])   && Corr_fact_fl == 'F') { iss >> dummy >> dummy >> conn_demand; }
		if (find_str("SEG LENGTH", Corr_fact_data[i]) && Corr_fact_fl == 'S') { iss >> dummy >> dummy >> seg_length; }
		
		if (find_str("FLOW CORR", Corr_fact_data[i])) {
			iss >> dummy >> dummy >> flow_corr_fl; if (flow_corr_fl == 'Y') { net->DE_options.flow_corr_fl = 1; }
		}
		if (find_str("DISP CORR", Corr_fact_data[i])) {
			iss >> dummy >> dummy >> disp_corr_fl; if (disp_corr_fl == 'Y') { net->DE_options.disp_corr_fl = 1; }
		}
		if (find_str("Rw CORR", Corr_fact_data[i])) {
			iss >> dummy >> dummy >> Rw_corr_fl; if (Rw_corr_fl == 'Y') { net->DE_options.Rw_corr_fl = 1; }
		}

	}

	if (Corr_fact_fl == 'F') {
		// display output
		WRITE_LOG_MSG("	o	Correction factors turned ON with a connection demand of " + toString(conn_demand) + " " + toString(net->options.Flow_UNITS));

		// change flow unit to m3/sec
		conn_demand *= net->options.Flow_unit_conv;

		// save the data
		net->DE_options.conn_demand = conn_demand;
	}
	else if (Corr_fact_fl == 'S') {
		// display output
		string length_units;

		if (net->options.unit_sys == 1) { length_units = "m"; } else { length_units = "ft"; }

		WRITE_LOG_MSG("	o	Correction factors turned ON with a segment spacing of " + toString(seg_length) + " " + toString(length_units));

		// Change spacing units to meters
		if (net->options.unit_sys == 0) { seg_length *= 0.3048; }
		
		// Save the data
		net->DE_options.seg_length = seg_length;
	}
	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Read Stochastic demand parameters
    vector<string> Stoc_dem_data;

	Stoc_dem_data = InputData(index, WUDESIMinp, "[STOC_DEMANDS]");

	for (int i = 0; i < Stoc_dem_data.size(); i++) {

		istringstream iss(Stoc_dem_data[i]);
		string dummy;
		char Stoc_dem_fl='N';

		// Read stochastic demands flag
		if (find_str("STOCHASTIC DEMANDS", Stoc_dem_data[i])) {
			
			iss >> dummy >> dummy >> Stoc_dem_fl;
			
			if (Stoc_dem_fl=='Y') {
				net->DE_options.Stoc_dem_fl = 1;
			}
			else if (Stoc_dem_fl=='N') {
				net->DE_options.Stoc_dem_fl = 0;
			}
			else {
				WRITE_LOG_MSG("	o	STOCHASTIC DEMANDS flag not recognized ");  return 1;
			}
		}

		// Read stochastic demands parameters
		if (find_str("u1", Stoc_dem_data[i])) { iss >> dummy >> net->DE_options.u1; }
		if (find_str("u2", Stoc_dem_data[i])) { iss >> dummy >> net->DE_options.u2; }
		if (find_str("s1", Stoc_dem_data[i])) { iss >> dummy >> net->DE_options.s1; }
		if (find_str("s2", Stoc_dem_data[i])) { iss >> dummy >> net->DE_options.s2; }
		if (find_str("AVERAGING INTERVAL", Stoc_dem_data[i])) { iss >> dummy >> dummy >> net->DE_options.avg_int; }

	}
	
	double dt_h_inp = net->times.Hyd_step_hr  * 3600. + net->times.Hyd_step_min  * 60.;  //User Input Hydraulic time step(sec)
	double dt_q     = net->times.Qual_step_hr * 3600. + net->times.Qual_step_min * 60.;  //User Input Quality time step(sec)

	if (net->DE_options.Stoc_dem_fl) {

		WRITE_LOG_MSG("	o	Stochastic demands turned ON with an averaging interval of " + toString(net->DE_options.avg_int) + " seconds");

		if (dt_h_inp <= net->DE_options.avg_int) { WRITE_LOG_MSG("Hydraulic step is not greater than the AVERAGING INTERVAL"); return 1; }
		if (dt_q     >= net->DE_options.avg_int) { WRITE_LOG_MSG("Water Quality step is not smaller than the AVERAGING INTERVAL"); return 1; }

	}
	else {
		WRITE_LOG_MSG("	o	Stochastic demands turned OFF ");
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read WUDESIM output options
	vector<string> WUDESIM_rpt_data;
	vector<double> simulated_branches;

	WUDESIM_rpt_data = InputData(index, WUDESIMinp, "[WUDESIM_OUTPUT]");

	char WUDESIM_rpt_fl = 'N';

	for (int i = 0; i < WUDESIM_rpt_data.size(); i++) {

		istringstream iss(WUDESIM_rpt_data[i]);
		string dummy;
		int branch_no = 0;

		// Read simulate all flag
		if (find_str("SIM ALL", WUDESIM_rpt_data[i])) {

			iss >> dummy >> dummy >> WUDESIM_rpt_fl;

			if (WUDESIM_rpt_fl == 'Y') {
				for (int j = 0; j < net->DE_branches.size(); ++j) { simulated_branches.push_back(j); }
			}
			else if (WUDESIM_rpt_fl != 'N') {
				WRITE_LOG_MSG("	o	SIM ALL flag not recognized!");  return 1;
			}
		}		

		// Read the IDs of the branches selected for simulation
		if (find_str("Branches", WUDESIM_rpt_data[i]) && WUDESIM_rpt_fl == 'N') {

			iss >> dummy;

			while (!iss.eof()) {
				iss >> branch_no;
				if (
					// check branch no is not zero
					(branch_no != 0) && 
					
					// check branch no is not out of the limit
					(branch_no <= net->DE_branches.size()) && 
					
					// check branch no is not already in the list
					(std::find(simulated_branches.begin(), simulated_branches.end(), branch_no-1) == simulated_branches.end())
					)
				{
					simulated_branches.push_back(branch_no-1); 
				}
			}
		}
	}	
	
	// store the data to the network
	net->DE_options.simulated_branches = simulated_branches;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;

}