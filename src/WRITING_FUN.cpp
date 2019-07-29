/*
Project:     WUDESIM ver. 1 BETA
File:        DEFIND.cpp
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
#include <iomanip>


#include "WUDESIM.h"
#include "Classes.h"
#include "WUDESIM_CORE.h"
#include "Utilities.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Define template printing function

template<typename T> void printElement(T t, ofstream& filename)
{
	filename << left << setw(14) << setfill(' ') << t;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write dead-end branch ids 
void write_DE_ids(Network* net)
{
	ofstream DE_IDS_FILE;

	if (!DE_IDS_FILE.is_open()) { DE_IDS_FILE.open("DE_Pipe_ID.out", ios::out | ios::trunc); }

	printElement("Branch_No", DE_IDS_FILE);
	printElement("Pipe_ID", DE_IDS_FILE);


	DE_IDS_FILE << endl;

	for (int branch = 0; branch < net->DE_branches.size(); branch++) {

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {
			
			printElement(branch + 1, DE_IDS_FILE);
			printElement(net->DE_branches[branch].pipe_id[DeadEnd], DE_IDS_FILE);
			DE_IDS_FILE << endl;

		}
	}
	DE_IDS_FILE.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_DE_Properties(Network* net) {

	ofstream DE_PROPERTIES_FILE;

	if (!DE_PROPERTIES_FILE.is_open()) { DE_PROPERTIES_FILE.open("DE_Properties.out", ios::out | ios::trunc); }

	printElement("Branch_No", DE_PROPERTIES_FILE);
	printElement("Pipe_ID", DE_PROPERTIES_FILE);
	printElement("Length(m)", DE_PROPERTIES_FILE);
	printElement("Diameter(m)", DE_PROPERTIES_FILE);
	printElement("Avg_Reynolds", DE_PROPERTIES_FILE);
	printElement("Avg_Res_T", DE_PROPERTIES_FILE);
	printElement("Avg_C_in", DE_PROPERTIES_FILE);
	printElement("Avg_C_out", DE_PROPERTIES_FILE);

	DE_PROPERTIES_FILE << endl;

	for (int branch = 0; branch < net->DE_branches.size(); branch++) {

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

			printElement(branch + 1, DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].pipe_id[DeadEnd], DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].length[DeadEnd], DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].diameter[DeadEnd], DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].Reynolds[DeadEnd]), DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].Res_time[DeadEnd]), DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].boundary_C_EPANET[DeadEnd]), DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].terminal_C_EPANET[DeadEnd]), DE_PROPERTIES_FILE);


			DE_PROPERTIES_FILE << endl;

		}
	}

	DE_PROPERTIES_FILE.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write stochastically generated demands
void write_stoc_dems(Network* net)
{
	ofstream DE_STOC_DEM_FILE;

	if (!DE_STOC_DEM_FILE.is_open()) { DE_STOC_DEM_FILE.open("DE_Stochastic_Flow.out", ios::out | ios::trunc);}

	// Find which branches will be simulated
	vector<double> DE_branch_simulation = net->DE_options.simulated_branches;

	for (int bb = 0; bb < DE_branch_simulation.size(); bb++) {

		int branch = DE_branch_simulation[bb];

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

			DE_STOC_DEM_FILE << "Stochastic Flow for Branch NO. " << branch + 1 << '\t' << "Pipe NO. " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;
			printElement("Timestep", DE_STOC_DEM_FILE);
			printElement("Flow (m3/s)", DE_STOC_DEM_FILE);
			DE_STOC_DEM_FILE << endl;

			for (int r_step = 0; r_step < net->DE_branches[branch].pipe_flow_WUDESIM[DeadEnd].size(); ++r_step) {
				printElement(r_step, DE_STOC_DEM_FILE);
				printElement(net->DE_branches[branch].pipe_flow_WUDESIM[DeadEnd][r_step], DE_STOC_DEM_FILE);
				DE_STOC_DEM_FILE << endl;
			}

			DE_STOC_DEM_FILE << endl;
		}
	}

	DE_STOC_DEM_FILE.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write WUDESIM report file
void write_WUDESIM_rpt(Network* net) {

	ofstream DE_REPORT_FILE;

	string WUDESIMrptfileNanme = net->WUDESIM_RPT;

	if (!DE_REPORT_FILE.is_open()) { DE_REPORT_FILE.open(WUDESIMrptfileNanme, ios::out | ios::trunc); }

	// Find which branches will be simulated
	vector<double> DE_branch_simulation = net->DE_options.simulated_branches;


	DE_REPORT_FILE << "*********************************************************************************************" << endl;
	DE_REPORT_FILE << "*                                                                                           *" << endl;
	DE_REPORT_FILE << "*                                     WUDESIM REPORT                                        *" << endl;
	DE_REPORT_FILE << "*                                                                                           *" << endl;
	DE_REPORT_FILE << "*********************************************************************************************" << endl << endl;

	// Read EPANET times variables
	double N_steps_EPANET = net->times.N_steps;                                     //Number of report steps
	double N_steps_skip = net->times.N_steps - net->times.N_steps_rep;            //Number of steps to skip when writing the report
	double dt_h_EPANET = net->times.Hyd_step_hr + net->times.Hyd_step_min / 60.; //Hydraulic time step(hr)

	// Set actual times variables
	double dt_h_act = net->DE_options.avg_int / 3600.;  // Averaging interval (sec->hr)
	int N_avg_int = dt_h_EPANET / dt_h_act;           // Number of averaging intervals per EPANET hydraulic step
	int N_steps_act = N_steps_EPANET * N_avg_int;       // Number of simulation hydraulic steps		

	for (int bb = 0; bb < DE_branch_simulation.size(); bb++) {

		int branch = DE_branch_simulation[bb];

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Get WUDESIM concentrations to write

			vector<double> WUDESIM_conc_write(N_steps_EPANET, 0);
			vector<double> WUDESIM_Pe_write(N_steps_EPANET, 0);

			if (net->DE_options.Stoc_dem_fl) {

				//Stochastic demand time steps
				int kk = 0;

				// Step through EPANET steps
				for (int epanet_step = 0; epanet_step < N_steps_EPANET; ++epanet_step) {

					// Store WUDESIM concentration directly at the end of each EPANET step
					WUDESIM_conc_write[epanet_step] = net->DE_branches[branch].terminal_C_WUDESIM[DeadEnd][kk];

					// Create dummy peclet number to calculate the average over EPANET steps
					vector<double> Peclet_dummy(N_avg_int, 0.);

					// Step through stochastic demand steps within each EPANET step
					for (int interv = 0; interv < N_avg_int; ++interv) {

						// Store the dummy Peclet number
						Peclet_dummy[interv] = net->DE_branches[branch].Peclet[DeadEnd][kk];

						++kk;
					}

					// Calculate the average peclet number within the EPANET step
					WUDESIM_Pe_write[epanet_step] = avrg(Peclet_dummy);
				}
			}
			else {
				WUDESIM_conc_write = net->DE_branches[branch].terminal_C_WUDESIM[DeadEnd];
				WUDESIM_Pe_write = net->DE_branches[branch].Peclet[DeadEnd];
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			

			// Write the output file
			printElement("Branch_No", DE_REPORT_FILE);
			printElement("Pipe_ID", DE_REPORT_FILE);
			printElement("Junc_ID", DE_REPORT_FILE);
			printElement("N_segments", DE_REPORT_FILE);
			printElement("Flow_Corr", DE_REPORT_FILE);
			printElement("Disp_Corr", DE_REPORT_FILE);
			printElement("Rw_Corr", DE_REPORT_FILE);
			printElement("Avg_Peclet", DE_REPORT_FILE);
			DE_REPORT_FILE << endl;

			printElement(branch + 1, DE_REPORT_FILE);
			printElement(net->DE_branches[branch].pipe_id[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].terminal_id[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].N_segment[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].Flow_Correction_factor[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].Disp_Correction_factor[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].Rw_Correction_factor[DeadEnd], DE_REPORT_FILE);
			printElement(avrg(net->DE_branches[branch].Peclet[DeadEnd]), DE_REPORT_FILE);
			DE_REPORT_FILE << endl;
			DE_REPORT_FILE << endl;

			DE_REPORT_FILE << "Simulated " << net->options.QUAL_TAG << " concnetrations (" << net->options.QUAL_UNIT << ")" << endl;
			printElement("Time", DE_REPORT_FILE);
			printElement("EPANET", DE_REPORT_FILE);
			printElement("WUDESIM", DE_REPORT_FILE);
			printElement("Reynolds", DE_REPORT_FILE);
			printElement("Peclet", DE_REPORT_FILE);
			printElement("Res_T", DE_REPORT_FILE);
			DE_REPORT_FILE << endl;

			for (int j = 0; j < N_steps_EPANET; ++j) {
				if (j >= N_steps_skip) {
					printElement(j * dt_h_EPANET, DE_REPORT_FILE);
					printElement(net->DE_branches[branch].terminal_C_EPANET[DeadEnd][j], DE_REPORT_FILE);
					printElement(WUDESIM_conc_write[j], DE_REPORT_FILE);
					printElement(net->DE_branches[branch].Reynolds[DeadEnd][j], DE_REPORT_FILE);
					printElement(WUDESIM_Pe_write[j], DE_REPORT_FILE);
					printElement(net->DE_branches[branch].Res_time[DeadEnd][j], DE_REPORT_FILE);
					DE_REPORT_FILE << endl;
				}
			}
			DE_REPORT_FILE << endl;
			DE_REPORT_FILE << "*********************************************************************************************" << endl << endl;
		}
	}

	DE_REPORT_FILE.close();
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////