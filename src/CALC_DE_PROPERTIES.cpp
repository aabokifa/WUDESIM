/*
Project:     WUDESIM ver. 1 BETA
File:        WQSIM.cpp
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


int CALC_DE_PROPERTIES(Network* net) {

	// Read EPANET time variables
	int N_steps_EPANET = net->times.N_steps;                                                //Number of hydraulic steps
	double dt_h_EPANET = net->times.Hyd_step_hr + net->times.Hyd_step_min / 60.;            //Hydraulic time step(hr)
	double dt_q        = net->times.Qual_step_hr * 3600. + net->times.Qual_step_min * 60.;  //User Input Quality time step(sec)

	// Set reactions variables
	double Kb  = net->reactions.Bulk_coeff;     //first order bulk decay coefficient(/sec)
	double Kw  = net->reactions.Wall_coeff;     //first order wall decay coefficeint(m/sec) or (1/m2/sec)
	double n_b = net->reactions.Bulk_order;	    //Bulk reaction order
	double n_w = net->reactions.Wall_order;		//Wall reaction order
	double C_L = net->reactions.Lim_pot;		//Limiting potential

	// Set transport variables
	double D_diff    = net->options.Rel_Diffusivity * 1.2E-9;      //molecular diffusion coefficient(m2 / sec)
	double viscosity = net->options.Rel_Viscosity * 1E-6;          //water kinematic viscosity(m2 / sec) = 1cSt

	for (int branch = 0; branch < net->DE_branches.size(); branch++) {

		int N_pipes = net->DE_branches[branch].branch_size;

		// Initialize property vectors
		net->DE_branches[branch].Reynolds.resize(N_pipes, vector<double>(N_steps_EPANET, 0.));
		net->DE_branches[branch].Res_time.resize(N_pipes, vector<double>(N_steps_EPANET, 0.));

		for (int DeadEnd = (N_pipes - 1); DeadEnd >= 0; DeadEnd--) {

			// read pipe flow
			vector<double> pipe_flow = net->DE_branches[branch].pipe_flow_EPANET[DeadEnd];
			
			// Pipe data
			double dp = net->DE_branches[branch].diameter[DeadEnd]; //pipe diameter (m)
			double Lt = net->DE_branches[branch].length[DeadEnd];   //pipe length (m)
			double r0 = dp / 2.;                                    //pipe radius(m)
			double Ap = 3.14159 * pow(r0, 2.);                      //pipe x - sec area(m2)

			// Flow velocity
			vector<double> u(N_steps_EPANET);
			vector<double> u_corr(N_steps_EPANET);
			
			for (int i = 0; i < N_steps_EPANET; ++i) {
				u[i]      = pipe_flow[i] / Ap;   //Actual Flow velocity(m / sec)
			}

			// Calculate Reynolds number
			for (int i = 0; i < N_steps_EPANET; ++i) {
				net->DE_branches[branch].Reynolds[DeadEnd][i] = u[i] * dp / viscosity;
			}

			// Calculate Residence Time
			for (int i = 0; i < N_steps_EPANET; ++i) {
				net->DE_branches[branch].Res_time[DeadEnd][i] = Lt / u[i];
			}
		}
	}


	return 0;
}