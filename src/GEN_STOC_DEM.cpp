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

int GEN_STOC_DEM(Network* net, vector<double> DE_branch_simulation)
{

	// EPANET demand parameters
	int N_steps_EPANET = net->times.N_steps;                                               //Number of hydraulic steps
	double dt_h_EPANET = net->times.Hyd_step_hr * 3600. + net->times.Hyd_step_min * 60.;   //EPANET Hydraulic time step (sec)
	
	// Stochastic pulse variables
	double u1      = net->DE_options.u1;	  //mean pulse duration[ln(s)]
	double s1      = net->DE_options.s1;	  //stdev pulse duration[ln(s)]
	double u2      = net->DE_options.u2;	  //mean pulse intensity[ln(L / s)]
	double s2      = net->DE_options.s2;	  //stdev pulse instenisty[ln(L / s)]
	double mean_V  = exp(u1 + u2) / 1000;	  // Mean pulse volume(m3)

	// Calculate new flow profile parameters
	double dt_h_WUDESIM             = net->DE_options.avg_int;		// Averaging interval (sec)
	int N_avg_int                   = dt_h_EPANET / dt_h_WUDESIM;   // Number of averaging intervals per EPANET hydraulic step
	int N_steps_WUDESIM             = N_steps_EPANET * N_avg_int;   // Number of WUDESIM hydraulic steps	
	net->DE_options.N_steps_WUDESIM = N_steps_WUDESIM;              // Store information to the network

	for (int bb = 0; bb < DE_branch_simulation.size(); bb++) {

		int branch = DE_branch_simulation[bb];

		int n_rows    = net->DE_branches[branch].branch_size;
		int n_columns = N_steps_WUDESIM;
		net->DE_branches[branch].pipe_flow_WUDESIM.resize(n_rows, vector<double>(n_columns, 0.));

		for (int DeadEnd = (net->DE_branches[branch].branch_size - 1); DeadEnd >= 0; DeadEnd--) {

			vector<double> flow_inp = net->DE_branches[branch].pipe_flow_EPANET[DeadEnd];
			
			// Calculate pulse arrival rate (Lambda)
			vector<double> L = flow_inp;
			for (int i = 0; i < N_steps_EPANET; ++i) { L[i] /= mean_V; }

			// Calculate number of pulses for each step
			vector<int> N_pulses(N_steps_EPANET, 0);
			int N_pulses_total = 0;
			for (int i = 0; i < N_steps_EPANET; ++i) { N_pulses[i] = round(dt_h_EPANET * L[i]); N_pulses_total += N_pulses[i]; }

			// Generate pulse durations and intensities
			int N_rows    = N_steps_EPANET;
			int N_columns = *max_element(N_pulses.begin(), N_pulses.end());

			vector<vector<double>> D, I;
			D.resize(N_rows, vector<double>(N_columns, 0.));
			I.resize(N_rows, vector<double>(N_columns, 0.));

			random_device seed_generator;
			default_random_engine generator(seed_generator());

			normal_distribution<double> D_distribution(u1, s1);
			normal_distribution<double> I_distribution(u2, s2);

			for (int step = 0; step < N_steps_EPANET; step++) {
				for (int pulse = 0; pulse < N_pulses[step]; pulse++) {
					D[step][pulse] = ceil(exp(D_distribution(generator)));
					I[step][pulse] = exp(I_distribution(generator));
				}
			}

			// Correct for total hourly demand volume
			vector<double> corr(N_steps_EPANET);
			for (int step = 0; step < N_steps_EPANET; step++) {
				double vol_sum = 0.;
				for (int pulse = 0; pulse < N_pulses[step]; pulse++) {
					vol_sum += (D[step][pulse] * I[step][pulse]);
				}
				corr[step] = (vol_sum / 1000) / (flow_inp[step] * dt_h_EPANET);
				for (int pulse = 0; pulse < N_pulses[step]; pulse++) { I[step][pulse] /= corr[step]; }
			}

			// Generate arrival times for pulses in each step
			vector<vector<double>> U, t;
			U.resize(N_rows, vector<double>(N_columns, 0.));
			t.resize(N_rows, vector<double>(N_columns, 0.));

			vector<double> Start(N_pulses_total);
			vector<double> End(N_pulses_total);
			vector<double> Intensity(N_pulses_total);
			vector<double> Duration(N_pulses_total);

			int k = 0;
			uniform_real_distribution<double> U_distribution(0.0, 1.0);
			for (int step = 0; step < N_steps_EPANET; step++) {
				for (int pulse = 0; pulse < N_pulses[step]; pulse++) {

					U[step][pulse] = U_distribution(generator);
					t[step][pulse] = ((-log(U[step][pulse])) / L[step]);

					Start[k] = accumulate(t[step].begin(), t[step].end(), 0) + dt_h_EPANET * step ;
					if (Start[k] < 0) { Start[k] == 0; }

					End[k] = Start[k] + D[step][pulse];
					if (End[k] > N_steps_EPANET * dt_h_EPANET ) { End[k] = N_steps_EPANET * dt_h_EPANET ; }

					Intensity[k] = I[step][pulse];
					Duration[k] = D[step][pulse];
					k++;
				}
			}

			// Sum flow from overlapping pulses
			vector<double> Q(N_steps_EPANET * dt_h_EPANET, 0); // Instantanous flows
			for (int pulse = 0; pulse < Start.size(); pulse++) {
				for (int i = Start[pulse]; i < End[pulse]; ++i) {
					Q[i] += Intensity[pulse];
				}
			}

			// Averaging flow rates
			vector<double> flow_act(N_steps_EPANET* N_avg_int, 0); // Averaged flows
			int T = 0, K = 0;
			for (int step = 0; step < N_steps_EPANET; step++) {
				for (int interv = 0; interv < N_avg_int; interv++) {
					for (int i = T; i < T + dt_h_WUDESIM; i++) { flow_act[K] += Q[i]; }
					flow_act[K] /= (dt_h_WUDESIM * 1000);
					T += dt_h_WUDESIM;
					K++;
				}
			}

			// Store averaged flows
			net->DE_branches[branch].pipe_flow_WUDESIM[DeadEnd] = flow_act;
		}
	}
	return 0;
}