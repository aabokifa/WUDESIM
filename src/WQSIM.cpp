/*
Project:     WUDESIM ver. 1 BETA
File:        WQSIM.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: This source file includes the water quality simulation function. It first reads the simulation inputs from the "WUDESIM.inp"
             file regarding the use of correction factors and stochastic flow demands. If the user choses to apply stochastic demands, 
		     the demand generation function, which resides in the "DEMGEN.cpp" source file, is called before the beginning of the water
			 quality simulation and the generated flow demands are written to the "WUDESIM.flow" file. The code then solves the 
			 advection-dispersion-reaction (ADR) equation for all dead-end branches in the network and the output is written to a new report 
			 file "WUDESIM.rpt", which includes the time series concentrations at the terminal junction of all dead-ends as simulated by
			 both EPANET and WUDESIM.
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


#include "Classes.h"
#include "WUDESIMmain.h"
#include "Utilities.h"

using namespace std;

ofstream myrptfile;
ofstream DEflow;

int WQSIM(string WUDESIMinpfileName, Network* net)
{
	// Import WUDESIM.inp file
	vector<string> WUDESIMinp;
	WUDESIMinp = ImportFile(WUDESIMinpfileName);

	if (WUDESIMinp.empty()) {
		cout << WUDESIMinpfileName << "is empty/corrupt!" << endl;
		return 1;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Find section locations
	vector<int> index;
	vector<string> Headers1 = { "[CORR_FACTS]","[STOC_DEMANDS]","[END]" };

	for (int j = 0;j < Headers1.size();++j) {
		for (int i = 0;i < WUDESIMinp.size();++i) { if (find(Headers1[j], WUDESIMinp[i])) { index.push_back(i); } };
	}
	sort(index.begin(), index.end());

	if (index.size() < 3) { cout << WUDESIMinpfileName << "must have [CORR_FACTS], [STOC_DEMANDS], and [END]" << endl; return 1; }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Read segment spacing
	double seg_spacing;
	bool Corr_Fact_fl;


	vector<string> Corr_fact_data;
	string dummy;
	Corr_fact_data = InputData(index, WUDESIMinp, "[CORR_FACTS]");

	for (int i = 0;i < Corr_fact_data.size(); ++i) {

		istringstream iss(Corr_fact_data[i]);

		if (find("CORRECTION FACTORS", Corr_fact_data[i])) {
			iss >> dummy >> dummy >> dummy;
			if (find("Y", dummy)) { Corr_Fact_fl = 1; }
			else { Corr_Fact_fl = 0; }
		}
		if (find("SEGMENT SPACING", Corr_fact_data[i])) {
			iss >> dummy >> dummy >> seg_spacing;
		}
	}	
	if (Corr_Fact_fl) {
		cout << "Correction Factors turned on with a segment spacing of " << seg_spacing << " meters" << endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Set times variables
	double N_steps_inp = net->times.N_steps;                                          //Number of report steps
	double dt_h_inp = net->times.Rep_step_hr + net->times.Rep_step_min / 60.;         //Hydraulic time step(hr)
	double dt_q = net->times.Qual_step_hr*3600. + net->times.Qual_step_min*60.;             //User Input Quality time step(sec)

	// Set reactions variables
	double Kb = net->reactions.Bulk_coeff;               //first order bulk decay coefficient(/sec)
	double Kw = net->reactions.Wall_coeff;               //first order wall decay coefficeint(m/sec) or (1/m2/sec)
	double n_b = net->reactions.Bulk_order;				 //Bulk reaction order
	double n_w = net->reactions.Wall_order;				 //Wall reaction order
	double C_L = net->reactions.Lim_pot;				 //Limiting potential
	
    // Set transport variables
	double D_diff = net->options.Rel_Diffusivity*1.2E-9;           //molecular diffusion coefficient(m2 / sec)
	double viscosity = net->options.Rel_Viscosity*1E-6;          //water kinematic viscosity(m2 / sec) = 1cSt

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Stochastic demand parameters
	bool Stoc_dem_fl;
	bool Taylor_fl;

	vector<double> Stoc_dem_pars(5);

	vector<string> Stoc_dem_data;

	Stoc_dem_data = InputData(index, WUDESIMinp, "[STOC_DEMANDS]");

	for (int i = 0;i < Stoc_dem_data.size();i++) {

		istringstream iss(Stoc_dem_data[i]);

		if (find("STOCHASTIC DEMANDS", Stoc_dem_data[i])) {
			iss >> dummy >> dummy >> dummy;
			if (find("Y", dummy)) {
				Stoc_dem_fl = 1;
				Taylor_fl = 0;
			}
			else {
				Stoc_dem_fl = 0;
				Taylor_fl = 1;
			}
		}

		if (find("u1", Stoc_dem_data[i])) { iss >> dummy >> Stoc_dem_pars[0]; }
		if (find("u2", Stoc_dem_data[i])) { iss >> dummy >> Stoc_dem_pars[1]; }
		if (find("s1", Stoc_dem_data[i])) { iss >> dummy >> Stoc_dem_pars[2]; }
		if (find("s2", Stoc_dem_data[i])) { iss >> dummy >> Stoc_dem_pars[3]; }
		if (find("AVERAGING INTERVAL", Stoc_dem_data[i])) { iss >> dummy >> dummy >> Stoc_dem_pars[4]; }

	}
	if (Stoc_dem_fl) {
		cout << "Stochastic demands turned on with an averaging interval of " << Stoc_dem_pars[4] << " seconds" << endl;
		if (dt_q >= Stoc_dem_pars[4]) { cout << "Water Quality step is not smaller than the Averaging Interval" << endl; return 1; }
		if (dt_h_inp * 3600 <= Stoc_dem_pars[4]) { cout << "Hydraulic step is not greater than the Averaging Interval" << endl; return 1; }
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Start Water Quality Simulations
	
	myrptfile.open("WUDESIM.rpt", ios::out | ios::trunc);
	DEflow.open("WUDESIM.flow", ios::out | ios::trunc);

	for (int branch = 0;branch < net->DE_branches.size();branch++) {

		int N_pipes = net->DE_branches[branch].branch_size;

		net->DE_branches[branch].terminal_new.resize(N_pipes, vector<double>(N_steps_inp, 0.));

		for (int DeadEnd = (N_pipes - 1); DeadEnd >= 0;DeadEnd--) {


			// Read Dead End Flow and boundary concentration
			vector<double> DEboundary(N_steps_inp, 0);
			vector<double> flow_inp(N_steps_inp,0);
				for (int r_step = 0;r_step < N_steps_inp;++r_step) {
					if (DeadEnd == N_pipes - 1) {
						DEboundary[r_step] = net->DE_branches[branch].boundary[DeadEnd][r_step];
					}
					else {
						DEboundary[r_step] = net->DE_branches[branch].terminal_new[DeadEnd + 1][r_step];
					}
					flow_inp[r_step]= net->DE_branches[branch].pipe_flow[DeadEnd][r_step];
				}			

			vector<double> flow_act;
			int N_steps_act;
			double dt_h_act;

			if (Stoc_dem_fl) {
				flow_act = DEMGEN(flow_inp, dt_q, dt_h_inp, Stoc_dem_pars);
				N_steps_act = flow_act.size();
				dt_h_act = N_steps_inp *dt_h_inp / N_steps_act;
			}
			else {
				flow_act = flow_inp;
				N_steps_act = N_steps_inp;
				dt_h_act = dt_h_inp;
			}

			double InitialConc = net->DE_branches[branch].boundary[DeadEnd][0];      // Initial quality is equal to upstream node
			vector<double> C_init;
			vector<double> C;
			vector<double> C_adv;
			vector<double> H;
			vector<double> GR;
			vector<double> GF;

			//Calculate Correction Factor for flow demand
			int N_seg_corr = 1;
			
			if (Corr_Fact_fl && DeadEnd == 0) {
				N_seg_corr = ceil(net->DE_branches[branch].length[DeadEnd] / seg_spacing);
			}

			double Flow_Corr = 0., Disp_Corr = 0., Kw_Corr = 0.;
			for (double i = 1;i < N_seg_corr + 1; ++i) {

				Flow_Corr = Flow_Corr + pow(N_seg_corr - i + 1., -1.);
				Disp_Corr = Disp_Corr + pow(N_seg_corr - i + 1., 2.);
				Kw_Corr = Kw_Corr + pow(N_seg_corr - i + 1., -2. / 3.);

			}

			Flow_Corr = 1. / Flow_Corr;
			Disp_Corr = Disp_Corr / (pow(N_seg_corr, 3.));
			Kw_Corr = Kw_Corr*Flow_Corr;


			// Pipe data
			double dp = net->DE_branches[branch].diameter[DeadEnd]; //pipe diameter (m)
			double Lt = net->DE_branches[branch].length[DeadEnd];   //pipe length (m)
			double r0 = dp / 2.;                  //pipe radius(m)
			double Ap = 3.14159*pow(r0, 2.);      //pipe x - sec area(m2)

			// Flow velocity
			vector<double> u(N_steps_act);
			vector<double> u_dum(N_steps_act);
			for (int i = 0;i < N_steps_act;++i) { u[i] = flow_act[i] / Ap; u_dum[i] = u[i] * Flow_Corr; };   //Actual and corrected Flow velocity(m / sec)

			//Check Quality time step is sufficient for max flow event
			double u_max = *max_element(u_dum.begin(), u_dum.end());     //maximum flow velocity(m / sec
			double dx_max = u_max*dt_q;                //maximum delta x(m)
			double N_min = Lt / dx_max + 1;            //min number of descritization points
			if (N_min < 10) {						   //Should have at least three descritization points
				N_min = 10;
				dx_max = Lt / (N_min - 1);
				dt_q = dx_max / u_max;
			}

			//Time discretization
			double Tot_time = dt_h_act*N_steps_act * 3600;        //Total time(sec)
			double Nqsteps = ceil(dt_h_act * 3600 / dt_q);        //No.of quality steps
			dt_q = dt_h_act * 3600 / Nqsteps;                     //Actual quality time step(sec)

			// Space discretaztion
			vector<double> N(N_steps_act);             //No. of descritization points in each hydraulic step
			for (int i = 0;i < N_steps_act;++i) { N[i] = floor(Lt / (u_dum[i] * dt_q)) + 1; if (isinf(N[i])) { N[i] = 1000; }; }
			vector<double> dx(N_steps_act);            //delta x of each hydraulic step
			for (int i = 0;i < N_steps_act;++i) { dx[i] = Lt / (N[i] - 1); }

			//Initialize terminal concentration
			net->DE_branches[branch].terminal_new[DeadEnd][0] = net->DE_branches[branch].terminal[DeadEnd][0];

			//Dynamic Dispersion Coefficient
			vector<double> E_taylor(N_steps_act);
			for (int i = 0;i < N_steps_act;++i) { E_taylor[i] = (pow(r0*u[i], 2) / (48. * D_diff))*Disp_Corr; }

			double** E_disp = new double*[Nqsteps];                                     //initialize the array
			for (int i = 0; i < Nqsteps; ++i) { E_disp[i] = new double[N_steps_act]; }  //initialize the array

			double T0 = pow(r0, 2.) / (16 * D_diff); //Eularian Time-scale

			for (int Hstep = 0; Hstep < N_steps_act;++Hstep) {     //now fill the array
				//Determine last nonzero flow event
				int i = Hstep;
				int k = 1;
				double E_prev = 0;
				double t_prev = 0;
				double u_prev = 0;

				while (k == 1) {
					if (i == 0) { k = 0; }
					else if (u[i - 1] > 0) {
						E_prev = E_disp[(int)Nqsteps - 1][i - 1];
						t_prev = (i*Nqsteps*dt_q);
						u_prev = u[i - 1];
						k = 0;
					}
					else { --i; k = 1; }
				}
				for (int qstep = 0; qstep < Nqsteps;++qstep) {
					double t_curr = (Hstep*Nqsteps + qstep + 1)*dt_q;
					if (Taylor_fl == 0) {
						if (u_prev == 0) { //First pulse
							E_disp[qstep][Hstep] = E_taylor[Hstep] * (1 - exp(-(qstep + 1)*dt_q / T0)) + D_diff;
						}
						else {             //Subsequent pulses
							E_disp[qstep][Hstep] = E_prev*(u[Hstep] / u_prev)*exp(-(t_curr - t_prev) / T0) + E_taylor[Hstep] * (1 - exp(-(qstep + 1)*dt_q / T0)) + D_diff;
						}
					}
					else if (Taylor_fl == 1) {
						E_disp[qstep][Hstep] = E_taylor[Hstep] + D_diff;
					}
				}
			}

			//Calculate dimensionless numbers and decay coefficient
			vector<double> Re(N_steps_act);
			vector<double> Pe(N_steps_act);
			vector<double> Sh(N_steps_act);
			vector<double> Kf(N_steps_act);
			vector<double> K(N_steps_act);
			double Sc = viscosity / D_diff; //Schmidt Number

			for (int i = 0;i < N_steps_act;++i) {
				Re[i] = u[i] * dp / viscosity;    //Reynolds Number
				Pe[i] = u[i] * Lt / E_taylor[i];  //Peclet Number based on Taylor's dispersion
				if (Re[i] >= 2300) { Sh[i] = 0.023*pow(Re[i], 0.83)*pow(Sc, 0.333); }   //Sherwood Numbr
				else { Sh[i] = 3.65 + (0.0668*(dp*Re[i] * Sc / Lt)) / (1 + 0.04*pow((dp*Re[i] * Sc / Lt), (2 / 3))); }
				Kf[i] = Sh[i] * D_diff / dp;                              //Mass transfer coefficient
			}

			// Boundary condition interpolation
			vector<double> C_bound_in(N_steps_inp + 1);
			vector<double> Times_bound_in(N_steps_inp + 1);

			for (int i = 0;i < N_steps_inp;++i) {
				C_bound_in[i] = DEboundary[i];
				Times_bound_in[i] = i*dt_h_inp;
			}
			C_bound_in[N_steps_inp] = C_bound_in[N_steps_inp - 1];    //Add another concentration to interpolate the last time step
			Times_bound_in[N_steps_inp] = N_steps_inp*dt_h_inp;


			vector<double> Times_bound_act(Nqsteps*N_steps_act);
			vector<double> C_bound_act(Nqsteps*N_steps_act);

			for (int i = 0;i < Times_bound_act.size();++i) { Times_bound_act[i] = i*dt_q / 3600; }
			C_bound_act = interpolation(Times_bound_in, C_bound_in, Times_bound_act);

			/////////////////////////////////////////////Loop for hydraulic time steps /////////////////////////////////////////

			int t = 0;

			for (int Hstep = 0;Hstep < N_steps_act;++Hstep) {

				//Space descritization
				vector<double> x(N[Hstep]);
				vector<double> E(N[Hstep]);

				for (int i = 0;i < N[Hstep];++i) {
					x[i] = i*dx[Hstep];        //Descritized space co - ordinate
					E[i] = x[i] - u_dum[Hstep] * dt_q; //Characteristic line footprint
				}
				double Q = flow_act[Hstep];          //flow rate of the current hydraulic step(m3 / sec)

				//Initial Condition
				if (t == 0) {
					C_init.resize(N[Hstep]);
					fill(C_init.begin(), C_init.end(), InitialConc);
				}
				C.clear();
				C.resize(N[Hstep]);
				C = C_init;

				/////////////////////////////////////////////Loop for quality time steps /////////////////////////////////////////

				for (int qstep = 0;qstep < Nqsteps;++qstep) {
					//Calculate Advection Concentration
					C_adv.clear();
					C_adv.resize(N[Hstep]);

					//Boundary Condition @x=0
					C[0] = C_bound_act[t];
					C_adv[0] = C_bound_act[t];

					if (u[Hstep] == 0.) { C_adv = C; }
					else {
						for (int i = 1;i < N[Hstep]; ++i) {
							//MOC advection step
							C_adv[i] = (C[i - 1] * (x[i] - E[i]) + C[i] * (E[i] - x[i - 1])) / dx[Hstep];

							//Bulk Reaction step
							if (C_L == 0) {
								if (n_b == 1) { 
									C_adv[i] *= exp(Kb* dt_q); }                           //First order bulk reaction
								else { C_adv[i] *= (1 + (n_b - 1)*Kb*pow(C_adv[i], n_b - 1)*dt_q); }   //nth order bulk reaction
							}
							else if (n_b > 0) { // dC/dt=Kb*(C_L-C)*C^(n_b-1)
								//Solve via RK4
								double y1 = C_adv[i];
								double k1 = Kb*(C_L - y1)*pow(y1, n_b - 1);
								double y2 = y1 + dt_q*k1 / 2;
								double k2 = Kb*(C_L - y2)*pow(y2, n_b - 1);
								double y3 = y1 + dt_q*k2 / 2;
								double k3 = Kb*(C_L - y3)*pow(y3, n_b - 1);
								double y4 = y1 + dt_q*k3;
								double k4 = Kb*(C_L - y4)*pow(y4, n_b - 1);

								C_adv[i] += (dt_q / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
							}
							else if (n_b < 0) {
								if (Kb < 0) {      // Michaelis-Menton Decay Kinetics: dC/dt=Kb*C/(C_L-C)
									double y1 = C_adv[i];
									double k1 = Kb*y1 / (C_L - y1);

									double y2 = y1 + dt_q*k1 / 2;
									double k2 = Kb*y2 / (C_L - y2);

									double y3 = y1 + dt_q*k2 / 2;
									double k3 = Kb*y3 / (C_L - y3);

									double y4 = y1 + dt_q*k3;
									double k4 = Kb*y4 / (C_L - y4);

									C_adv[i] += (dt_q / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

								}
								else if (Kb > 0) { //Michaelis-Menton growth Kinetics dC/dt=Kb*C/(C_L+C)
									double y1 = C_adv[i];
									double k1 = Kb*y1 / (C_L + y1);

									double y2 = y1 + dt_q*k1 / 2;
									double k2 = Kb*y2 / (C_L + y2);

									double y3 = y1 + dt_q*k2 / 2;
									double k3 = Kb*y3 / (C_L + y3);

									double y4 = y1 + dt_q*k3;
									double k4 = Kb*y4 / (C_L + y4);

									C_adv[i] += (dt_q / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
								}
							}
							//Wall Reaction 
							double Rw;
							if (n_w == 1) {        //First order wall reaction
								Rw = ((4. * Kw*Kf[Hstep] )/ (dp*(Kw + Kf[Hstep])))*Kw_Corr;
								if (Kw + Kf[Hstep] <= 0) { Rw = 0; }
								C_adv[i] *= exp(Rw* dt_q);
							}
							else if (n_w == 0) {   //Zeroth order wall reaction
								Rw = (4 * Kw / dp) / 1000; // Kw/rh (mg/m3/sec --> mg/L/sec)
								C_adv[i] += Rw*dt_q;
							}
						}
					}


					//Calculate lambda
					double  L = E_disp[qstep][Hstep] * dt_q / pow(dx[Hstep], 2);

					//Create the tridiagonal matrix (A)
					int n = N[Hstep] - 2;
					vector<double> a1(n);
					vector<double> a2(n);
					vector<double> a3(n);
					fill(a1.begin() + 1, a1.end(), -L);
					fill(a2.begin(), a2.end(), 1 + 2 * L);
					fill(a3.begin(), a3.end() - 1, -L);

					// Solve for H
					vector<double> H_str(n);
					for (int i = 0; i < n; ++i) { H_str[i] = C_adv[i + 1]; }
					solveThomas(a1, a2, a3, H_str, n);

					H.clear();
					H.resize(N[Hstep]);
					for (int i = 1; i <= n; ++i) { H[i] = H_str[i - 1]; }


					// Solve for GR
					vector<double> GR_str(n);
					GR_str.front() = L;
					solveThomas(a1, a2, a3, GR_str, n);

					GR.clear();
					GR.resize(N[Hstep]);
					GR.front() = 1.;
					for (int i = 1; i <= n; ++i) { GR[i] = GR_str[i - 1]; }

					// Solve for GF
					vector<double> GF_str(n);
					GF_str.back() = L;
					solveThomas(a1, a2, a3, GF_str, n);

					GF.clear();
					GF.resize(N[Hstep]);
					GF.back() = 1.;
					for (int i = 1; i <= n; ++i) { GF[i] = GF_str[i - 1]; }

					// Calculate concentrations
					for (int i = 0; i < N[Hstep];++i) { C[i] = H[i] + GR[i] * C_adv[0] + GF[i] * C_adv[N[Hstep] - 1]; }

					t++;

				} //Loop for Quality steps

				// Store terminal concentration
				if (Hstep*dt_h_act / dt_h_inp - floor(Hstep*dt_h_act / dt_h_inp) == 0 || Hstep > 0) {
					int T = Hstep*dt_h_act / dt_h_inp;
					net->DE_branches[branch].terminal_new[DeadEnd][T] = C[N[Hstep] - 1];
				}


				//Now interpolate the concentration profile to fit the new hydraulic step
				if (Hstep < N_steps_act - 1) {
					C_init.clear();
					C_init.resize(N[Hstep + 1]);
					vector<double> x_next(N[Hstep + 1]);

					for (int j = 0;j < N[Hstep + 1];++j) {
						x_next[j] = j*dx[Hstep + 1];
					}
					C_init = interpolation(x, C, x_next);

				}

			} //Loop for Hydraulic Steps


			// Write terminal concentration and generated flow to output file
			if (myrptfile.is_open()) {
				myrptfile << "Terminal "<<  net->options.QUAL_TAG << " concnetration (" << net->options.QUAL_UNIT << ") for branch no. " << branch + 1 << endl;
				myrptfile << "Junction No. " << net->DE_branches[branch].terminal_id[DeadEnd] << endl;
				myrptfile << "Time" << '\t' << "EPANET" << '\t' << "WUDESIM" << endl;

				for (int j = 0;j < N_steps_inp;++j) {
					myrptfile << j*dt_h_inp << '\t' << net->DE_branches[branch].terminal[DeadEnd][j] << '\t' << net->DE_branches[branch].terminal_new[DeadEnd][j] << endl;
				}
				myrptfile << endl;
			}

			if (DEflow.is_open() && Stoc_dem_fl) {
				DEflow << "Generated flow for branch no. " << branch + 1 << endl;
				DEflow << "Pipe No. " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;
				DEflow << "Time" << '\t' << "Generated Demand (m3/s)" << endl;

				for (int j = 0;j < N_steps_act;++j) {
					DEflow << j << '\t' << flow_act[j] <<  endl;
				}
				DEflow << endl;
			}

		}


	}
	return 0;
}

