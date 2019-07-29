/*
Project:     WUDESIM ver. 1 BETA
File:        DEMGEN.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
*/


#include <iostream> 
#include <fstream>  
#include <vector>   
#include <map>
#include <algorithm>
#include <string>
#include <random>
#include <numeric>

#include "WUDESIM.h"
#include "Classes.h"
#include "WUDESIM_CORE.h"
#include "Utilities.h"
#include "WRITING_FUN.h"

using namespace std;

int CALC_CORR_FACT(Network* net) {

	// Read correction factor parameters
	int    Corr_Fact_fl = net->DE_options.Corr_Fact_fl;
	double conn_demand  = net->DE_options.conn_demand;
	double seg_spacing  = net->DE_options.seg_length;

	for (int branch = 0; branch < net->DE_branches.size(); branch++) {

		int N_pipes = net->DE_branches[branch].branch_size;

		// Initialize correction factors
		net->DE_branches[branch].N_segment.resize(N_pipes, 1);
		net->DE_branches[branch].Flow_Correction_factor.resize(N_pipes, 1.);
		net->DE_branches[branch].Disp_Correction_factor.resize(N_pipes, 1.);
		net->DE_branches[branch].Rw_Correction_factor.resize(N_pipes, 1.);

		for (int DeadEnd = (N_pipes - 1); DeadEnd >= 0; DeadEnd--) {

			// Initialize Correction Factors
			double Flow_Corr = 1., Disp_Corr = 1., Kw_Corr = 1.;
			int    N_seg_corr = 1, N_seg_dum = 1;

			if (Corr_Fact_fl != 0) {

				// Compare pipe inflow and outflow
				vector<double> Q_in  = net->DE_branches[branch].pipe_flow_EPANET[DeadEnd];				
				double Qin_avg  = avrg(Q_in);
				double Qout_avg = 0.;
				
				if (DeadEnd != 0) { //For the last dead-end pipe: Qout_avg = 0
					vector<double> Q_out = net->DE_branches[branch].pipe_flow_EPANET[DeadEnd - 1];
					Qout_avg = avrg(Q_out);
				}
				
				// Calculate dummy number of segments
				if (Corr_Fact_fl == 1) {      //flow-based correction
					N_seg_dum = ceil((Qin_avg - Qout_avg) / conn_demand);
				}
				else if (Corr_Fact_fl == 2) { //spacing-based correction
					N_seg_dum = ceil(net->DE_branches[branch].length[DeadEnd] / seg_spacing);
					conn_demand = (Qin_avg - Qout_avg) / N_seg_dum;
				}

				// Apply the correction factors only if a significant drop in the flow takes place
				if (Qout_avg <= conn_demand) {

					Flow_Corr = 0.;
					Disp_Corr = 0.;
					Kw_Corr   = 0.;

					for (double i = 1; i <= N_seg_dum ; ++i) {
						Flow_Corr = Flow_Corr + pow(N_seg_dum - i + 1., -1.);
						Disp_Corr = Disp_Corr + pow(N_seg_dum - i + 1., 2.);
						Kw_Corr   = Kw_Corr   + pow(N_seg_dum - i + 1., -(2. / 3.));
					}

					N_seg_corr = N_seg_dum;
					Flow_Corr = 1. / Flow_Corr;
					Disp_Corr = Disp_Corr / (pow(N_seg_corr, 3.));
					Kw_Corr = Kw_Corr * Flow_Corr;
				}
				else {
					N_seg_corr = 1.;
					Flow_Corr = 1.;
					Disp_Corr = 1.;
					Kw_Corr = 1.;
				}
			}

			net->DE_branches[branch].N_segment[DeadEnd] = N_seg_corr;
			if (net->DE_options.flow_corr_fl) { net->DE_branches[branch].Flow_Correction_factor[DeadEnd] = Flow_Corr; }
			if (net->DE_options.disp_corr_fl) {	net->DE_branches[branch].Disp_Correction_factor[DeadEnd] = Disp_Corr; }
			if (net->DE_options.Rw_corr_fl) { net->DE_branches[branch].Rw_Correction_factor[DeadEnd]     = Kw_Corr; }
		}
	}

	return 0;
}