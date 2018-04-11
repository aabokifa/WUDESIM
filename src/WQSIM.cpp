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

int WQSIM(string WUDESIMinpfileName, string WUDESIMrptfileNanme, Network* net)
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
		cout << "o	Correction Factors turned on with a segment spacing of " << seg_spacing << " meters" << endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Set times variables
	double N_steps_inp = net->times.N_steps;                                          //Number of report steps
	double dt_h_inp = net->times.Hyd_step_hr + net->times.Hyd_step_min / 60.;         //Hydraulic time step(hr)
	double dt_q = net->times.Qual_step_hr*3600. + net->times.Qual_step_min*60.;       //User Input Quality time step(sec)
	double N_steps_skip = net->times.N_steps - net->times.N_steps_rep;				  //Number of steps to skip when writing the report

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
		cout << "o	Stochastic demands turned on with an averaging interval of " << Stoc_dem_pars[4] << " seconds" << endl;
		if (dt_q >= Stoc_dem_pars[4]) { cout << "Water Quality step is not smaller than the Averaging Interval" << endl; return 1; }
		if (dt_h_inp * 3600 <= Stoc_dem_pars[4]) { cout << "Hydraulic step is not greater than the Averaging Interval" << endl; return 1; }
		DEflow.open("WUDESIM.flow", ios::out | ios::trunc);
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Start Water Quality Simulations
	
	myrptfile.open(WUDESIMrptfileNanme, ios::out | ios::trunc);

	for (int branch = 0;branch < net->DE_branches.size();branch++) {

		cout << "o	Dead-End Branch no." << branch + 1 << endl;

		int N_pipes = net->DE_branches[branch].branch_size;

		net->DE_branches[branch].terminal_new.resize(N_pipes, vector<double>(N_steps_inp, 0.));
		net->DE_branches[branch].Reynolds.resize(N_pipes, vector<double>(N_steps_inp, 0.));
		net->DE_branches[branch].Correction_factors.resize(N_pipes, vector<double>(3, 1.));
		net->DE_branches[branch].N_segment.resize(N_pipes);

		for (int DeadEnd = (N_pipes - 1); DeadEnd >= 0;DeadEnd--) {

			dt_q = net->times.Qual_step_hr*3600. + net->times.Qual_step_min*60.;       //Reset WQ to the User Input Quality time step(sec)

			cout << "		Pipe ID " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;

			//Initialize terminal concentration
			net->DE_branches[branch].terminal_new[DeadEnd][0] = net->DE_branches[branch].terminal[DeadEnd][0];


			// Read Dead End Flow and boundary concentration
			vector<double> DEboundary(N_steps_inp, 0);
			vector<double> flow_inp(N_steps_inp, 0);
			for (int r_step = 0;r_step < N_steps_inp;++r_step) {
				if (DeadEnd == N_pipes - 1) {
					DEboundary[r_step] = net->DE_branches[branch].boundary[DeadEnd][r_step];
				}
				else {
					DEboundary[r_step] = net->DE_branches[branch].terminal_new[DeadEnd + 1][r_step];
				}
				flow_inp[r_step] = net->DE_branches[branch].pipe_flow[DeadEnd][r_step];
			}


			vector<double> flow_act;
			int N_steps_act;
			double dt_h_act;

			if (Stoc_dem_fl) {
				flow_act    = DEMGEN(flow_inp, dt_q, dt_h_inp, Stoc_dem_pars);
				N_steps_act = flow_act.size();
				dt_h_act    = N_steps_inp * dt_h_inp / N_steps_act;
			}
			else {
				flow_act = flow_inp;
				N_steps_act = N_steps_inp;
				dt_h_act = dt_h_inp;
			}

			double InitialConc = net->DE_branches[branch].terminal[DeadEnd][0];      // Initial quality is equal to downstream node
			vector<double> C_init;
			vector<double> C;
			vector<double> C_adv;
			vector<double> H;
			vector<double> GR;
			vector<double> GF;

			//Calculate Correction Factors
			double Flow_Corr = 1., Disp_Corr = 1., Kw_Corr = 1.;
			net->DE_branches[branch].N_segment[DeadEnd] = 1;

			if (Corr_Fact_fl) {

				int N_seg_corr = ceil(net->DE_branches[branch].length[DeadEnd] / seg_spacing);

				// Compare pipe inflow and outflow
				double Qin_avg = avrg(flow_inp);
				double Qout_avg = 0.;

				vector<double> flow_out(N_steps_inp, 0);

				if (DeadEnd != 0) { //For the last dead-end pipe: Qout_avg = 0

					for (int r_step = 0;r_step < N_steps_inp;++r_step) {
						flow_out[r_step] = net->DE_branches[branch].pipe_flow[DeadEnd - 1][r_step];
					}
					Qout_avg = avrg(flow_out);
				}

				if (Qout_avg <= (3 * Qin_avg / N_seg_corr)) {
					Flow_Corr = 0.;
					Disp_Corr = 0.;
					Kw_Corr   = 0.;
					
					for (double i = 1;i < N_seg_corr + 1; ++i) {
						Flow_Corr = Flow_Corr + pow(N_seg_corr - i + 1., -1.);
						Disp_Corr = Disp_Corr + pow(N_seg_corr - i + 1., 2.);
						Kw_Corr   = Kw_Corr   + pow(N_seg_corr - i + 1., -(2. / 3.));
					}

					Flow_Corr = 1. / Flow_Corr;
					Disp_Corr = Disp_Corr / (pow(N_seg_corr, 3.));
					Kw_Corr   = Kw_Corr * Flow_Corr;

					net->DE_branches[branch].N_segment[DeadEnd] = N_seg_corr;
					net->DE_branches[branch].Correction_factors[DeadEnd][0] = Flow_Corr;
					net->DE_branches[branch].Correction_factors[DeadEnd][1] = Disp_Corr;
					net->DE_branches[branch].Correction_factors[DeadEnd][2] = Kw_Corr;
				}
			}


			// Pipe data
			double dp = net->DE_branches[branch].diameter[DeadEnd]; //pipe diameter (m)
			double Lt = net->DE_branches[branch].length[DeadEnd];   //pipe length (m)
			double r0 = dp / 2.;                                    //pipe radius(m)
			double Ap = 3.14159*pow(r0, 2.);                        //pipe x - sec area(m2)

			// Flow velocity
			vector<double> u(N_steps_act);
			vector<double> u_dum(N_steps_act);
			for (int i = 0;i < N_steps_act;++i) {
				u[i] = flow_act[i] / Ap;     //Actual Flow velocity(m / sec)
				u_dum[i] = u[i] * Flow_Corr; //Corrected Flow velocity(m / sec)

				if (Lt / u_dum[i] > 3600 * 24 * 7.) { u_dum[i] = 0; } //Consider stagnant if residence time is more than 24 hrs


			}

			// Calculate Reynolds number
			vector<double> Re(N_steps_act);
			for (int i = 0;i < N_steps_act;++i) {
				Re[i] = u[i] * dp / viscosity;       //Calculate Reynolds Number
			}


			//Check Quality time step is sufficient for max flow event
			double u_max = abs(*max_element(u_dum.begin(), u_dum.end()));     //maximum flow velocity(m / sec)
			double dx_max = u_max * dt_q;              //maximum delta x(m)
			double N_min = Lt / dx_max + 1;            //min number of descritization points
			if (N_min < 10) {						   //Should have at least ten descritization points
				N_min = 10;
				dx_max = Lt / (N_min - 1);
				dt_q = dx_max / u_max;
			}

			//Time discretization
			double Tot_time = dt_h_act * N_steps_act * 3600;      //Total time(sec)
			double Nqsteps = ceil(dt_h_act * 3600 / dt_q);        //No.of quality steps
			dt_q = dt_h_act * 3600 / Nqsteps;                     //Actual quality time step(sec)

			// Space discretaztion
			vector<double> N(N_steps_act);             //No. of descritization points in each hydraulic step
			for (int i = 0;i < N_steps_act;++i) { 
				
				N[i] = floor(Lt / (u_dum[i] * dt_q)) + 1; 
							
				if (isinf(N[i])) { N[i] = 100; }; 
			
			}
			vector<double> dx(N_steps_act);            //delta x of each hydraulic step
			for (int i = 0;i < N_steps_act;++i) { dx[i] = Lt / (N[i] - 1); }
						
			//Dynamic Dispersion Coefficient
			vector<double> E_taylor(N_steps_act);

			for (int i = 0;i < N_steps_act;++i) {						
					E_taylor[i] = (pow(r0*u[i], 2) / (48. * D_diff))*Disp_Corr;
			}

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
							E_disp[qstep][Hstep] = E_prev * (u[Hstep] / u_prev)*exp(-(t_curr - t_prev) / T0) + E_taylor[Hstep] * (1 - exp(-(qstep + 1)*dt_q / T0)) + D_diff;
						}
					}
					else if (Taylor_fl == 1) {
						if (Re[Hstep] < 2300) {        //Laminar flow --> use Taylor's dispersion
							E_disp[qstep][Hstep] = E_taylor[Hstep] + D_diff;
						}
						else if (Re[Hstep] < 4000) { // Transitional flow --> interpolate to get dispersion							
							E_disp[qstep][Hstep] = E_taylor[Hstep] * (1 - (Re[Hstep] - 2300) / 1700) + D_diff;
						}
					}
				}
			}

			//Calculate dimensionless numbers and decay coefficient
			vector<double> Pe(N_steps_act);
			vector<double> Sh(N_steps_act);
			vector<double> Kf(N_steps_act);
			vector<double> K(N_steps_act);
			double Sc = viscosity / D_diff; //Schmidt Number
			double Re_avg = 0; //Avg Reynolds number
			double Pe_avg = 0; //Avg Peclet number

			for (int i = 0;i < N_steps_act;++i) {

				//Calculate Average Reynolds number
				Re_avg = Re_avg + Re[i];
				net->DE_branches[branch].Reynolds[DeadEnd][i] = Re[i];
				

				//Peclet Number based on Taylor's dispersion
				Pe[i] = u[i] * Lt / E_taylor[i];
				Pe_avg = Pe_avg + Pe[i];

				//Sherwood Numbr
				if (Re[i] >= 2300) { 
				
					// For turbulent flow --> corrected Reynolds number is used to give reduced mass transfer
					Sh[i] = 0.023*pow(Re[i]*Flow_Corr, 0.83)*pow(Sc, 0.333); 
				
				}
				else { 
					
					// For laminar flow --> original Reynolds number is used because the wall demand will later be corrected
					if (Corr_Fact_fl) {
						// use segment length
						int N_seg_corr = ceil(net->DE_branches[branch].length[DeadEnd] / seg_spacing);
						double Lseg = Lt / N_seg_corr;
						Sh[i] = 3.65 + (0.0668*(dp*Re[i] * Sc / Lseg)) / (1 + 0.04*pow((dp*Re[i] * Sc / Lseg), (2 / 3)));
					}
					else {
						// use full pipe length
						Sh[i] = 3.65 + (0.0668*(dp*Re[i] * Sc / Lt)) / (1 + 0.04*pow((dp*Re[i] * Sc / Lt), (2 / 3)));
					}
				}

				//Mass transfer coefficient
				Kf[i] = Sh[i] * D_diff / dp;
			}
			Re_avg = Re_avg / N_steps_act;
			Pe_avg = Pe_avg / N_steps_act;

			// Boundary condition interpolation
			vector<double> C_bound_in(N_steps_inp + 1);
			vector<double> Times_bound_in(N_steps_inp + 1);

			for (int i = 0;i < N_steps_inp;++i) {
				C_bound_in[i] = DEboundary[i];
				Times_bound_in[i] = i * dt_h_inp;
			}
			C_bound_in[N_steps_inp] = C_bound_in[N_steps_inp - 1];    //Add another concentration to interpolate the last time step
			Times_bound_in[N_steps_inp] = N_steps_inp * dt_h_inp;


			vector<double> Times_bound_act(Nqsteps*N_steps_act);
			vector<double> C_bound_act(Nqsteps*N_steps_act);

			for (int i = 0;i < Times_bound_act.size();++i) { Times_bound_act[i] = i * dt_q / 3600; }
			C_bound_act = interpolation(Times_bound_in, C_bound_in, Times_bound_act);

			/////////////////////////////////////////////Loop for hydraulic time steps /////////////////////////////////////////

			int t = 0;

			for (int Hstep = 0;Hstep < N_steps_act - 1;++Hstep) {


				//Space descritization
				vector<double> x(N[Hstep]);
				vector<double> E(N[Hstep]);

				for (int i = 0;i < N[Hstep];++i) {
					x[i] = i * dx[Hstep];                 //Descritized space co - ordinate
					E[i] = x[i] - u_dum[Hstep] * dt_q;    //Characteristic line footprint
				}
				double Q = flow_act[Hstep];             //flow rate of the current hydraulic step(m3 / sec)


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


					// Advection Step
					C_adv.clear();
					C_adv.resize(N[Hstep]);

					//Boundary Condition @x=0
					C[0] = C_bound_act[t];
					C_adv[0] = C_bound_act[t];

					if (u[Hstep] == 0.) {
						C_adv = C;
					}
					else {
						for (int i = 1;i < N[Hstep]; ++i) {
							//MOC advection step
							C_adv[i] = (C[i - 1] * (x[i] - E[i]) + C[i] * (E[i] - x[i - 1])) / dx[Hstep];


							//Bulk Reaction step
							if (C_L == 0) {
								if (n_b == 1) {
									C_adv[i] *= exp(Kb* dt_q);
								}                                       //First order bulk reaction
								else { C_adv[i] *= (1 + (n_b - 1)*Kb*pow(C_adv[i], n_b - 1)*dt_q); }   //nth order bulk reaction
							}
							else if (n_b > 0) { // dC/dt=Kb*(C_L-C)*C^(n_b-1)
								//Solve via RK4
								double y1 = C_adv[i];
								double k1 = Kb * (C_L - y1)*pow(y1, n_b - 1);
								double y2 = y1 + dt_q * k1 / 2;
								double k2 = Kb * (C_L - y2)*pow(y2, n_b - 1);
								double y3 = y1 + dt_q * k2 / 2;
								double k3 = Kb * (C_L - y3)*pow(y3, n_b - 1);
								double y4 = y1 + dt_q * k3;
								double k4 = Kb * (C_L - y4)*pow(y4, n_b - 1);

								C_adv[i] += (dt_q / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
							}
							else if (n_b < 0) {
								if (Kb < 0) {      // Michaelis-Menton Decay Kinetics: dC/dt=Kb*C/(C_L-C)
									double y1 = C_adv[i];
									double k1 = Kb * y1 / (C_L - y1);

									double y2 = y1 + dt_q * k1 / 2;
									double k2 = Kb * y2 / (C_L - y2);

									double y3 = y1 + dt_q * k2 / 2;
									double k3 = Kb * y3 / (C_L - y3);

									double y4 = y1 + dt_q * k3;
									double k4 = Kb * y4 / (C_L - y4);

									C_adv[i] += (dt_q / 6)*(k1 + 2 * k2 + 2 * k3 + k4);

								}
								else if (Kb > 0) { //Michaelis-Menton growth Kinetics dC/dt=Kb*C/(C_L+C)
									double y1 = C_adv[i];
									double k1 = Kb * y1 / (C_L + y1);

									double y2 = y1 + dt_q * k1 / 2;
									double k2 = Kb * y2 / (C_L + y2);

									double y3 = y1 + dt_q * k2 / 2;
									double k3 = Kb * y3 / (C_L + y3);

									double y4 = y1 + dt_q * k3;
									double k4 = Kb * y4 / (C_L + y4);

									C_adv[i] += (dt_q / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
								}
							}
							//Wall Reaction 
							double Rw;
							if (n_w == 1) {        //First order wall reaction
								
								if (Re[Hstep] < 2300) {
									// include correction factor only for laminar flow
									Rw = ((4. * abs(Kw)*Kf[Hstep]) / (dp*(abs(Kw) + Kf[Hstep])))*Kw_Corr;
								}
								else
								{
									Rw = ((4. * abs(Kw)*Kf[Hstep]) / (dp*(abs(Kw) + Kf[Hstep])));
								}


								if (Kw < 0) {
									C_adv[i] *= exp(-Rw* dt_q);
								}
								else if (Kw > 0){
									C_adv[i] *= exp(Rw* dt_q);
								}								
							}
							else if (n_w == 0) {            // Zeroth order wall reaction
								Rw = (4 * Kw / dp) / 1000;  // Kw/rh (mg/m3/sec --> mg/L/sec)
								C_adv[i] += Rw * dt_q;
							}
						}
					}

					// Dispersion Step (Only for Laminar & Transitional flows)
					if (Re[Hstep] < 4000)
					{
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
					}
					else
					{
						// No dispersion for turbulent flow
						for (int i = 0; i < N[Hstep];++i) {
							C[i] = C_adv[i];
						}
					}

					t++;


				} //Loop for Quality steps


				// Store terminal concentration
				if (Hstep*dt_h_act / dt_h_inp - floor(Hstep*dt_h_act / dt_h_inp) == 0) {
					int T = Hstep * dt_h_act / dt_h_inp + 1;
					net->DE_branches[branch].terminal_new[DeadEnd][T] = C[N[Hstep] - 1];
					
					if (abs(net->DE_branches[branch].terminal_new[DeadEnd][T])<0.0001) {
						net->DE_branches[branch].terminal_new[DeadEnd][T] = 0;
					}
				}



				//Now interpolate the concentration profile to fit the new hydraulic step
				if (Hstep < N_steps_act - 1) {
					C_init.clear();
					C_init.resize(N[Hstep + 1]);
					vector<double> x_next(N[Hstep + 1]);

					for (int j = 0;j < N[Hstep + 1];++j) {
						x_next[j] = j * dx[Hstep + 1];
					}
					C_init = interpolation(x, C, x_next);
				}



			} //Loop for Hydraulic Steps


			// Write terminal concentration and generated flow to output file
			if (myrptfile.is_open()) {
				myrptfile << "Terminal " << net->options.QUAL_TAG << " concnetration (" << net->options.QUAL_UNIT << ") for branch no. " << branch + 1 << endl;
				myrptfile << "Pipe ID: " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;
				myrptfile << "Junction ID: " << net->DE_branches[branch].terminal_id[DeadEnd] << endl;
				myrptfile << "Avg Re number = " << Re_avg << " ; No. of segments = " << net->DE_branches[branch].N_segment[DeadEnd] <<" ; Length = "<< net->DE_branches[branch].length[DeadEnd] << " ; Diameter = " << net->DE_branches[branch].diameter[DeadEnd] << endl;
				myrptfile << "Flow_Corr_Fact = " << net->DE_branches[branch].Correction_factors[DeadEnd][0] << " ; Disp_Corr_Fact = " << net->DE_branches[branch].Correction_factors[DeadEnd][1] << " ; Rw_Corr_Fact = " << net->DE_branches[branch].Correction_factors[DeadEnd][2] << endl;
				myrptfile << "Time" << '\t' << "EPANET" << '\t' << "WUDESIM" << '\t' << "Reynolds" << endl;

				for (int j = 0;j < N_steps_inp;++j) {
					if (j >= N_steps_skip) {
						myrptfile
							<< j * dt_h_inp
							<< '\t' << net->DE_branches[branch].terminal[DeadEnd][j]
							<< '\t' << net->DE_branches[branch].terminal_new[DeadEnd][j]
							<< '\t' << net->DE_branches[branch].Reynolds[DeadEnd][j]
							<< endl;
					}
				}
				myrptfile << endl;
			}

			if (DEflow.is_open() && Stoc_dem_fl) {
				DEflow << "Generated flow for branch no. " << branch + 1 << endl;
				DEflow << "Pipe No. " << net->DE_branches[branch].pipe_id[DeadEnd] << endl;
				DEflow << "Time" << '\t' << "Generated Demand (m3/s)" << endl;

				for (int j = 0;j < N_steps_act;++j) {
					DEflow << j << '\t' << flow_act[j] << endl;
				}
				DEflow << endl;
			}



		}
	}
	return 0;
}

