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
#include "WRITING_FUN.h"
#include "Utilities.h"
#include "epanet2.h"


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

	if (!DE_IDS_FILE.is_open()) { DE_IDS_FILE.open("DEB_ID.OUT", ios::out | ios::trunc); }

	DE_IDS_FILE << "*********************************************************************************************" << endl;
	DE_IDS_FILE << "*                                                                                           *" << endl;
	DE_IDS_FILE << "*                              IDs of Dead End Branch Pipes                                 *" << endl;
	DE_IDS_FILE << "*                                                                                           *" << endl;
	DE_IDS_FILE << "*********************************************************************************************" << endl << endl;

	printElement("Branch_ID", DE_IDS_FILE);
	printElement("Branch_Size", DE_IDS_FILE);
	printElement("Pipe_IDs", DE_IDS_FILE);

	DE_IDS_FILE << endl;

	for (int branch = 0; branch < net->DE_branches.size(); branch++) {

		printElement(net->DE_branches[branch].branch_id, DE_IDS_FILE);
		printElement(net->DE_branches[branch].branch_size, DE_IDS_FILE);

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {			
			printElement(net->DE_branches[branch].pipe_id[DeadEnd], DE_IDS_FILE);
		}

		DE_IDS_FILE << endl;

	}

	DE_IDS_FILE << "*********************************************************************************************" << endl << endl;

	DE_IDS_FILE.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_EPANET_rpt(Network* net) {
	ENreport();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_DE_Properties(Network* net) {

	ofstream DE_PROPERTIES_FILE;

	if (!DE_PROPERTIES_FILE.is_open()) { DE_PROPERTIES_FILE.open("DEB_PROP.OUT", ios::out | ios::trunc); }

	DE_PROPERTIES_FILE << "*********************************************************************************************" << endl;
	DE_PROPERTIES_FILE << "*                                                                                           *" << endl;
	DE_PROPERTIES_FILE << "*                            Properties of Dead End Branches                                *" << endl;
	DE_PROPERTIES_FILE << "*                                 from EPANET Simulation                                    *" << endl;
	DE_PROPERTIES_FILE << "*                                                                                           *" << endl;
	DE_PROPERTIES_FILE << "*********************************************************************************************" << endl << endl;

	printElement("Branch_ID", DE_PROPERTIES_FILE);
	printElement("Pipe_ID", DE_PROPERTIES_FILE);
	printElement("Length(m)", DE_PROPERTIES_FILE);
	printElement("Diameter(m)", DE_PROPERTIES_FILE);
	printElement("Avg_Reynolds", DE_PROPERTIES_FILE);
	printElement("Avg_Res_T", DE_PROPERTIES_FILE);
	printElement("Junc_ID", DE_PROPERTIES_FILE);
	printElement("Avg_QUAL", DE_PROPERTIES_FILE);

	DE_PROPERTIES_FILE << endl;


	for (int branch = 0; branch < net->DE_branches.size(); branch++) {

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

			printElement(net->DE_branches[branch].branch_id, DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].pipe_id[DeadEnd], DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].length[DeadEnd], DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].diameter[DeadEnd], DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].Reynolds_EPANET[DeadEnd]), DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].Res_time_EPANET[DeadEnd]), DE_PROPERTIES_FILE);
			printElement(net->DE_branches[branch].terminal_id[DeadEnd], DE_PROPERTIES_FILE);
			printElement(avrg(net->DE_branches[branch].terminal_C_EPANET[DeadEnd]), DE_PROPERTIES_FILE);

			DE_PROPERTIES_FILE << endl;

		}
	}

	DE_PROPERTIES_FILE << "*********************************************************************************************" << endl << endl;

	DE_PROPERTIES_FILE.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Write stochastically generated demands
void write_stoc_dems(Network* net)
{
	ofstream DE_STOC_DEM_FILE;

	if (!DE_STOC_DEM_FILE.is_open()) { DE_STOC_DEM_FILE.open("DEB_STOC_FL.OUT", ios::out | ios::trunc);}

	DE_STOC_DEM_FILE << "*********************************************************************************************" << endl;
	DE_STOC_DEM_FILE << "*                                                                                           *" << endl;
	DE_STOC_DEM_FILE << "*                             Stochastically Generated Flows                                *" << endl;
	DE_STOC_DEM_FILE << "*                                                                                           *" << endl;
	DE_STOC_DEM_FILE << "*********************************************************************************************" << endl << endl;

	// Find which branches will be simulated
	vector<double> DE_branch_simulation = net->DE_options.simulated_branches;

	for (int bb = 0; bb < DE_branch_simulation.size(); bb++) {

		int branch = DE_branch_simulation[bb];

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

			DE_STOC_DEM_FILE << "Stochastic Flow for Branch " << net->DE_branches[branch].branch_id << '\t' << "Pipe " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;
			printElement("Timestep", DE_STOC_DEM_FILE);
			printElement("Flow (m3/s)", DE_STOC_DEM_FILE);
			DE_STOC_DEM_FILE << endl;

			for (int r_step = 0; r_step < net->DE_branches[branch].pipe_flow_STOC[DeadEnd].size(); ++r_step) {
				printElement(r_step, DE_STOC_DEM_FILE);
				printElement(net->DE_branches[branch].pipe_flow_STOC[DeadEnd][r_step], DE_STOC_DEM_FILE);
				DE_STOC_DEM_FILE << endl;
			}

			DE_STOC_DEM_FILE << endl;
		}
	}

	DE_STOC_DEM_FILE << "*********************************************************************************************" << endl << endl;

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


	for (int bb = 0; bb < DE_branch_simulation.size(); bb++) {

		int branch = DE_branch_simulation[bb];

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			

			// Write the output file
			printElement("Branch_ID", DE_REPORT_FILE);
			printElement("Pipe_ID", DE_REPORT_FILE);
			printElement("Junc_ID", DE_REPORT_FILE);
			printElement("N_segments", DE_REPORT_FILE);
			printElement("Flow_Corr", DE_REPORT_FILE);
			printElement("Disp_Corr", DE_REPORT_FILE);
			printElement("Rw_Corr", DE_REPORT_FILE);
			DE_REPORT_FILE << endl;

			printElement(net->DE_branches[branch].branch_id, DE_REPORT_FILE);
			printElement(net->DE_branches[branch].pipe_id[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].terminal_id[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].N_segment[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].Flow_Correction_factor[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].Disp_Correction_factor[DeadEnd], DE_REPORT_FILE);
			printElement(net->DE_branches[branch].Rw_Correction_factor[DeadEnd], DE_REPORT_FILE);
			DE_REPORT_FILE << endl;
			DE_REPORT_FILE << endl;

			DE_REPORT_FILE << '\t' << net->options.QUAL_TAG << "(" << net->options.QUAL_UNIT << ") Results";
			printElement("", DE_REPORT_FILE);
			printElement("", DE_REPORT_FILE);
			DE_REPORT_FILE << '\t' << "Pipe Results" << endl;

			printElement("Time", DE_REPORT_FILE);
			printElement("EPANET", DE_REPORT_FILE);
			printElement("WUDESIM", DE_REPORT_FILE);
			printElement("", DE_REPORT_FILE);
			printElement("Reyn_WUDESIM", DE_REPORT_FILE);
			printElement("Pecl_WUDESIM", DE_REPORT_FILE);
			printElement("ResT_WUDESIM", DE_REPORT_FILE);
			DE_REPORT_FILE << endl;

			for (int j = 0; j < N_steps_EPANET; ++j) {
				if (j >= N_steps_skip) {
					printElement(j * dt_h_EPANET, DE_REPORT_FILE);
					printElement(net->DE_branches[branch].terminal_C_EPANET[DeadEnd][j], DE_REPORT_FILE);
					printElement(net->DE_branches[branch].terminal_C_WUDESIM_WRITE[DeadEnd][j], DE_REPORT_FILE);
					printElement("", DE_REPORT_FILE);
					printElement(net->DE_branches[branch].Reynolds_WUDESIM_WRITE[DeadEnd][j], DE_REPORT_FILE);
					printElement(net->DE_branches[branch].Peclet_WUDESIM_WRITE[DeadEnd][j], DE_REPORT_FILE);
					printElement(net->DE_branches[branch].Res_time_WUDESIM_WRITE[DeadEnd][j], DE_REPORT_FILE);
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