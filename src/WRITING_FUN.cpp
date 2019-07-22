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
	ofstream ofs;
	ofs.open("DE_Pipe_ID.out", ios::out | ios::trunc);
	if (ofs.is_open()) {

		ofs << "Deadend pipe IDs:" << endl;

		for (int branch = 0; branch < net->DE_branches.size(); branch++) {
			ofs << "DE Branch No. " << branch + 1 << " :" << endl;

			for (int pipe = 0; pipe < net->DE_branches[branch].branch_size; pipe++) {
				ofs << net->DE_branches[branch].pipe_id[pipe] << endl;
			}
		}
		ofs.close();
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_DE_Properties(Network* net) {

	ofstream outfile;

	outfile.open("DE_Properties.out", ios::out | ios::trunc);

	if (outfile.is_open()) {

		printElement("Branch_No", outfile);
		printElement("Pipe_ID", outfile);
		printElement("Length(m)", outfile);
		printElement("Diameter(m)", outfile);
		printElement("Avg_Reynolds", outfile);
		printElement("Avg_Res_T", outfile);
		printElement("Avg_C_in", outfile);
		printElement("Avg_C_out", outfile);
		printElement("C_out/C_in", outfile);

		outfile << endl;

		for (int branch = 0; branch < net->DE_branches.size(); branch++) {

			for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

				printElement(branch + 1, outfile);
				printElement(net->DE_branches[branch].pipe_id[DeadEnd], outfile);
				printElement(net->DE_branches[branch].length[DeadEnd], outfile);
				printElement(net->DE_branches[branch].diameter[DeadEnd], outfile);
				printElement(avrg(net->DE_branches[branch].Reynolds[DeadEnd]), outfile);
				printElement(avrg(net->DE_branches[branch].Res_time[DeadEnd]), outfile);
				printElement(avrg(net->DE_branches[branch].boundary_C_EPANET[DeadEnd]), outfile);
				printElement(avrg(net->DE_branches[branch].terminal_C_EPANET[DeadEnd]), outfile);
				printElement(avrg(net->DE_branches[branch].terminal_C_EPANET[DeadEnd])/
					         avrg(net->DE_branches[branch].boundary_C_EPANET[DeadEnd]), outfile);

				outfile << endl;

			}
		}
		outfile << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write stochastically generated demands
void write_stoc_dems(Network* net)
{
	// Find which branches will be simulated
	vector<double> DE_branch_simulation = net->DE_options.simulated_branches;

	ofstream DEflow;

	DEflow.open("DE_Stochastic_Flow.out", ios::out | ios::trunc);

	if (DEflow.is_open()) {

		for (int bb = 0; bb < DE_branch_simulation.size(); bb++) {

			int branch = DE_branch_simulation[bb];

			for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

				DEflow << "Stochastic Flow for Branch NO. " << branch + 1 << '\t' << "Pipe NO. " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;
				printElement("Timestep", DEflow);
				printElement("Flow (m3/s)", DEflow);
				DEflow << endl;

				for (int r_step = 0; r_step < net->DE_branches[branch].pipe_flow_WUDESIM[DeadEnd].size(); ++r_step) {
					printElement(r_step, DEflow);
					printElement(net->DE_branches[branch].pipe_flow_WUDESIM[DeadEnd][r_step], DEflow);
					DEflow << endl;
				}

				DEflow << endl;
			}

		}
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write WUDESIM report file
void write_WUDESIM_rpt(Network* net) {
	
	// Find which branches will be simulated
	vector<double> DE_branch_simulation = net->DE_options.simulated_branches;

	string WUDESIMrptfileNanme = net->WUDESIM_RPT;

	ofstream myrptfile;

	myrptfile.open(WUDESIMrptfileNanme, ios::out | ios::trunc);

	if (myrptfile.is_open()) {
		
		myrptfile << "*********************************************************************************************" << endl;
		myrptfile << "*                                                                                           *" << endl;
		myrptfile << "*                                     WUDESIM REPORT                                        *" << endl;
		myrptfile << "*                                                                                           *" << endl;
		myrptfile << "*********************************************************************************************" << endl << endl;

		// Read EPANET times variables
		double N_steps_EPANET  = net->times.N_steps;                                     //Number of report steps
		double N_steps_skip    = net->times.N_steps - net->times.N_steps_rep;            //Number of steps to skip when writing the report
		double dt_h_EPANET     = net->times.Hyd_step_hr + net->times.Hyd_step_min / 60.; //Hydraulic time step(hr)
		
		// Set actual times variables
		double dt_h_act = net->DE_options.avg_int / 3600.;  // Averaging interval (sec->hr)
		int N_avg_int   = dt_h_EPANET / dt_h_act;           // Number of averaging intervals per EPANET hydraulic step
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
				printElement("Branch_No", myrptfile);
				printElement("Pipe_ID", myrptfile);
				printElement("Junc_ID", myrptfile);
				printElement("N_segments", myrptfile);
				printElement("Flow_Corr", myrptfile);
				printElement("Disp_Corr", myrptfile);
				printElement("Rw_Corr", myrptfile);
				printElement("Avg_Peclet", myrptfile);				
				myrptfile << endl;
				
				printElement(branch + 1, myrptfile);
				printElement(net->DE_branches[branch].pipe_id[DeadEnd], myrptfile);
				printElement(net->DE_branches[branch].terminal_id[DeadEnd], myrptfile);
				printElement(net->DE_branches[branch].N_segment[DeadEnd], myrptfile);
				printElement(net->DE_branches[branch].Correction_factors[DeadEnd][0], myrptfile);
				printElement(net->DE_branches[branch].Correction_factors[DeadEnd][1], myrptfile);
				printElement(net->DE_branches[branch].Correction_factors[DeadEnd][2], myrptfile);
				printElement(avrg(net->DE_branches[branch].Peclet[DeadEnd]), myrptfile);
				myrptfile << endl;				
				myrptfile << endl;

				myrptfile << "Simulated " << net->options.QUAL_TAG << " concnetrations (" << net->options.QUAL_UNIT << ")" << endl;
				printElement("Time", myrptfile);
				printElement("EPANET", myrptfile);
				printElement("WUDESIM", myrptfile);
				printElement("Reynolds", myrptfile);
				printElement("Peclet", myrptfile);
				printElement("Res_T", myrptfile);			
				myrptfile << endl;

				for (int j = 0; j < N_steps_EPANET; ++j) {
					if (j >= N_steps_skip) {
						printElement(j* dt_h_EPANET, myrptfile);
						printElement(net->DE_branches[branch].terminal_C_EPANET[DeadEnd][j], myrptfile);
						printElement(WUDESIM_conc_write[j], myrptfile);
						printElement(net->DE_branches[branch].Reynolds[DeadEnd][j], myrptfile);
						printElement(WUDESIM_Pe_write[j], myrptfile);
						printElement(net->DE_branches[branch].Res_time[DeadEnd][j], myrptfile);
						myrptfile << endl;
					}
				}
				myrptfile << endl;
				myrptfile << "*********************************************************************************************" << endl << endl;
			}
		}
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////