/*
File:   WQSIM.cpp
Author: Ahmed Abokifa
Date:   10/25/2016
Desc:	This source file includes the water quality simulation function, which reads its inputs (pipe lengths and diameters, flow rates, 
		boundary and terminal concentrations, and transport and reaction parameters) from the binary files written by "DEFIND.cpp". It 
		then reads the simulation inputs specified by the user in the "WUDESIM.inp" file regarding the use of correction factors and 
		stochastic flow demands. If the user choses to apply stochastic demands, the demand generation function, which resides in the 
		"DEMGEN.cpp" source file, is called before the beginning of the water quality simulation. The code then solves the 
		advection-dispersion-reaction (ADR) equation for all dead-ends in the network. The output is written to a separate report file 
		"WUDESIM.rpt", which includes the time series concentrations at the terminal junction of all dead-ends as simulated by both 
		EPANET and WUDESIM.
*/


//#include "stdafx.h" 
#include <iostream> 
#include <fstream>  
#include <vector>  
#include <map>
#include <algorithm>
#include <sstream>  
#include <string>  
#include <iterator>
#include <numeric>
#include "Headers.h"
#include <stdio.h>
#include <tchar.h>
#include <SDKDDKVer.h>

using namespace std;



int WUDESIM(int branch_no, vector<string>& term_id)
{

	ifstream  ifs;

	//Read Simulation Inputs

	vector<double> Inputs(10);
	ifs.open("SIMinfo.bin", ios::in | ios::binary);
	if (ifs.is_open()) {
		ifs.read(reinterpret_cast<char*>(Inputs.data()), 7 * (sizeof(double)));
		ifs.close();
	}
	else {
		logfile << "Error opening SIMinfo.bin file" << endl;
	}

	double N_steps_inp = Inputs[0];
	double dt_h_inp = Inputs[1];         //Hydraulic time step(hr)
	double dt_q = Inputs[2];             //User Input Quality time step(sec)
	double Kb = Inputs[3];               //first order bulk decay coefficient(/sec)
	double Kw = Inputs[4];               //first order wall decay coefficeint(m/sec) or (1/m2/sec)
	double D_diff = Inputs[5];           //molecular diffusion coefficient(m2 / sec)
	double viscosity = Inputs[6];        //water kinematic viscosity(m2 / sec) = 1cSt
	double n_b = Inputs[7];				 //Bulk reaction order
	double n_w = Inputs[8];				 //Wall reaction order
	double C_L = Inputs[9];				 //Limiting potential

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Read Dead End branch info

	int N_pipes;	//Number of pipes in dead end branch
	vector<double> DE_lengths;
	vector<double> DE_diameters;

	ifs.open("DEinfo.bin", ios::in | ios::binary);
	if (ifs.is_open()) {

		ifs.read(reinterpret_cast<char*>(&N_pipes), sizeof(int)); //First read number of pipes in DE branch

		DE_lengths.resize(N_pipes);
		DE_diameters.resize(N_pipes);

		for (int i = 0;i < N_pipes;++i) {
			ifs.read(reinterpret_cast<char*>(&DE_lengths[i]), sizeof(double));
			ifs.read(reinterpret_cast<char*>(&DE_diameters[i]), sizeof(double));
		}
		ifs.close();
	}	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read Dead End flows

	int n_rows = N_pipes;
	int n_columns = N_steps_inp;

	vector<vector<double>> DEflow;
	DEflow.resize(n_rows, vector<double>(n_columns, 0.));

	ifs.open("DEflow.bin", ios::in | ios::binary);
	if (ifs.is_open()) {
		for (int j = 0;j < N_pipes;++j) {
			for (int r_step = 0;r_step < N_steps_inp;++r_step) {
				ifs.read(reinterpret_cast<char*>(&DEflow[j][r_step]), sizeof(double));				
			}
		}
		ifs.close();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Read Dead End terminal concentrations

	vector<vector<double>> DEterminal;
	DEterminal.resize(n_rows, vector<double>(n_columns, 0.));

	ifs.open("DEterminal.bin", ios::in | ios::binary);
	if (ifs.is_open()) {
		for (int j = 0;j < N_pipes;++j) {
			for (int r_step = 0;r_step < N_steps_inp;++r_step) {
				ifs.read(reinterpret_cast<char*>(&DEterminal[j][r_step]), sizeof(double));
			}
		}
		ifs.close();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Read WUDESIM.inp file

	vector<string> WUDESIMinp;
	WUDESIMinp = ImportFile("WUDESIM.inp");
	vector<string> Headers1 = { "[INP_FILE]","[RPT_FILE]","[CORR_FACTS]","[STOC_DEMANDS]","[END]" };

	vector<int> index1;
	for (int j = 0;j < Headers1.size();++j) {
		for (int i = 0;i < WUDESIMinp.size();++i) { if (find(Headers1[j], WUDESIMinp[i])) { index1.push_back(i); } };
	}
	sort(index1.begin(), index1.end());

	vector<string> dum;

	vector<int> N_seg_corr(N_pipes,1.0);
	double seg_spacing;
	dum = InputData(index1, WUDESIMinp, "[CORR_FACTS]");
	if (find("Y", dum[0])) {
		istringstream iss(dum[1]);
		iss >> seg_spacing;
		N_seg_corr[0] = ceil(DE_lengths[0] / seg_spacing); //Correction factors only apply to the dead end link
		logfile << "Correction factors will be used in the simulations" << endl;
		logfile << "No. of corrected segments for branch " << branch_no + 1 << " is " << N_seg_corr[0] << "segments" << endl; 
	}


	double Stoc_dem_fl;
	double Taylor_fl;

	vector<double> Stoc_dem_pars(5);

	string dummy;
	dum = InputData(index1, WUDESIMinp, "[STOC_DEMANDS]");

	if (find("Y", dum[0])) {
		Stoc_dem_fl = 1;
		Taylor_fl = 0;

		for (int i = 1;i <= 5;i++) {

			istringstream iss(dum[i]);

			if (find("u1", dum[i])) { iss >> dummy >> Stoc_dem_pars[0]; }
			if (find("u2", dum[i])) { iss >> dummy >> Stoc_dem_pars[1]; }
			if (find("s1", dum[i])) { iss >> dummy >> Stoc_dem_pars[2]; }
			if (find("s2", dum[i])) { iss >> dummy >> Stoc_dem_pars[3]; }
			if (find("Avg_int", dum[i])) { iss >> dummy >> Stoc_dem_pars[4]; }

		}
	}
	else { Stoc_dem_fl = 0; Taylor_fl = 1; }


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		//Start Water Quality Simulations

	for (int DeadEnd = (N_pipes - 1); DeadEnd >= 0;DeadEnd--) {


		// Read Dead End boundary concentration
		vector<double> DEboundary(N_steps_inp, 0);

		ifs.open("DEboundary.bin", ios::in | ios::binary);
		if (ifs.is_open()) {
			for (int r_step = 0;r_step < N_steps_inp;++r_step) {
				ifs.read(reinterpret_cast<char*>(&DEboundary[r_step]), sizeof(double));

			}
			ifs.close();		
		}

	    
		vector<double> flow_inp(N_steps_inp);
		for (int i = 0;i < N_steps_inp;i++) { flow_inp[i] = DEflow[DeadEnd][i]; }

		vector<double> flow_act;
		int N_steps_act;
		double dt_h_act;

		if (Stoc_dem_fl == 1) {
			flow_act = DEMGEN(flow_inp, dt_q, dt_h_inp, Stoc_dem_pars);
			N_steps_act = flow_act.size();
			dt_h_act = N_steps_inp *dt_h_inp / N_steps_act;
		}
		else {
			flow_act = flow_inp;
			N_steps_act = N_steps_inp;
			dt_h_act = dt_h_inp;
		}

		double InitialConc = (DEboundary[0] + DEterminal[DeadEnd][0]) / 2;      //Initial quality(mg / L)
		vector<double> C_init;
		vector<double> C;
		vector<double> C_adv;
		vector<double> C_terminal;
		vector<double> H;
		vector<double> GR;
		vector<double> GF;

		//Calculate Correction Factor for flow demand
		double Flow_Corr = 0., Disp_Corr = 0., Kw_Corr = 0.;

		for (double i = 1;i < N_seg_corr[DeadEnd] + 1; ++i) {

			Flow_Corr = Flow_Corr + pow(N_seg_corr[DeadEnd] - i + 1., -1.);
			Disp_Corr = Disp_Corr + pow(N_seg_corr[DeadEnd] - i + 1., 2.);
			Kw_Corr = Kw_Corr + pow(N_seg_corr[DeadEnd] - i + 1., -2. / 3.);

		}

		Flow_Corr = 1. / Flow_Corr;
		Disp_Corr = Disp_Corr / (pow(N_seg_corr[DeadEnd], 3.));
		Kw_Corr = Kw_Corr*Flow_Corr;


		// Pipe data
		double dp = DE_diameters[DeadEnd];
		double Lt = DE_lengths[DeadEnd];
		double r0 = dp / 2.;             //pipe radius(m)
		double Ap = 3.14159*pow(r0, 2.);      //pipe x - sec area(m2)

		// Flow velocity
		vector<double> u(N_steps_act);
		vector<double> u_dum(N_steps_act);
		for (int i = 0;i < N_steps_act;++i) { u[i] = flow_act[i] / Ap; u_dum[i] = u[i] * Flow_Corr; };   //Actual and corrected Flow velocity(m / sec)

		//Check Quality time step is sufficient for max flow event
		double u_max = *max_element(u_dum.begin(), u_dum.end());     //maximum flow velocity(m / sec
		double dx_max = u_max*dt_q;                //maximum delta x(m)
		double N_min = Lt / dx_max + 1;            //min number of descritization points
		if (N_min < 3) {						   //Should have at least three descritization points
			N_min = 3;
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
		C_terminal.resize(N_steps_inp);
		C_terminal[0] = InitialConc;

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
			Pe[i] = u[i] * Lt / E_taylor[i]; //Peclet Number based on Taylor's dispersion
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
		C_bound_in[N_steps_inp] = C_bound_in[N_steps_inp - 1]; //Add another concentration to interpolate the last time step
		Times_bound_in[N_steps_inp] = N_steps_inp*dt_h_inp;


		vector<double> Times_bound_act(Nqsteps*N_steps_act);
		vector<double> C_bound_act(Nqsteps*N_steps_act);

		for (int i = 0;i < Times_bound_act.size();++i) { Times_bound_act[i] = i*dt_q / 3600; }
		C_bound_act = interpolation(Times_bound_in, C_bound_in, Times_bound_act);

		/////////////////////////////////////////////Loop for hydraulic time steps /////////////////////////////////////////

		int t = 0;
		int T = 0;

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
							if (n_b == 1) { C_adv[i] *= exp(Kb* dt_q); }                           //First order bulk reaction
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
							Rw = (4 * Kw*Kf[i] / (dp*(Kw + Kf[i])))*Kw_Corr;
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
			
			if (T / dt_h_act - floor(T / dt_h_act) == 0) {
				C_terminal[T]= C[N[Hstep] - 1];
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
			T++;

		} //Loop for Hydraulic Steps


		// Write terminal concentration to output file
		if (myrptfile.is_open()) {
			myrptfile << "Terminal concentrations for branch no. " << branch_no+1 << endl;
			myrptfile << "Junction No. " << term_id[DeadEnd] << endl;
			myrptfile << "Time" << '\t' << "EPANET" << '\t' << "WUDESIM" << endl;

			for (int j = 0;j < N_steps_inp;++j) {
				myrptfile << j*dt_h_inp << '\t' << DEterminal[DeadEnd][j] << '\t' << C_terminal[j] << endl;
			}
			myrptfile << endl;
		}
		else { logfile << "rpt file not written!" << endl; }

		//Write boundary condition for next pipe
		ofstream  ofs;
		ofs.open("DEboundary.bin", ios::out | ios::binary | ios::trunc);
		if (ofs.is_open()) {
			for (int step = 0;step < N_steps_inp;++step) {
				ofs.write(reinterpret_cast<char*>(&C_terminal[step]), sizeof(double));
			}

			ofs.close();
		}
	}

	logfile << "Simulation for branch no. " << branch_no+1 << " done!" << endl << endl;






	return 0;
}

