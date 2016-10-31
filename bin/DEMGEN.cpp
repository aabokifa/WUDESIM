/*
File:   WQSIM.cpp
Author: Ahmed Abokifa
Date:   10/25/2016
Desc:	This file contains the stochastic demand generation function, which can be called from the water quality simulator. This function takes 
		the generalized flow rates, as well as the probability distribution parameters specified by the user in "WUDESIM.inp" in order to generate 
		stochastic flow demands based on Poisson pulse distributions at fine time scales. The stochastic demands have randomly generated intensities 
		and durations based on lognormal probability distributions. The demands are then aggregated and averaged based on the averaging interval 
		specified by the user. The function returns a vector of the stochastically generated demands to be used in the water quality simulations. 
*/



//#include "stdafx.h"
#include <iostream> 
#include <fstream>  
#include <vector>   
#include <map>
#include <algorithm>
#include <string>
#include <random>
#include <numeric>
#include <stdio.h>
#include <tchar.h>
#include <SDKDDKVer.h>

using namespace std;

vector<double> DEMGEN(vector<double>& flow_inp, double dt_q, double dt_h_inp, vector<double>& Stoc_dem_pars)
{

	// Stochastic demand variables
	double Avg_int = Stoc_dem_pars[4];		//Averaging interval in seconds
	
	double u1 = Stoc_dem_pars[0];			//mean pulse duration[ln(s)]
	double s1 = Stoc_dem_pars[2];			//stdev pulse duration[ln(s)]
	
	double u2 = Stoc_dem_pars[1];			//mean pulse intensity[ln(L / s)]
	double s2 = Stoc_dem_pars[3];			//stdev pulse instenisty[ln(L / s)]
	
	double mean_V = exp(u1+u2)/1000;	//Mean pulse volume(m3)

	int N_steps_inp = flow_inp.size();

	if (dt_q >= Avg_int) { cout << "Water Quality step is not smaller than the Averaging Interval" << endl; }
	if (dt_h_inp*3600 <= Avg_int) { cout << "Hydraulic step is not greater than the Averaging Interval" << endl; }

	// Calculate pulse arrival rate (Lambda)
	vector<double> L = flow_inp;
	for (int i = 0;i < N_steps_inp;++i) { L[i] /= mean_V; }

	// Calculate number of pulses for each step
	vector<int> N_pulses(N_steps_inp,0);
	int N_pulses_total = 0;
	for (int i = 0;i < N_steps_inp;++i) { N_pulses[i] = round(dt_h_inp * 3600 * L[i]); N_pulses_total += N_pulses[i];}
	//cout << N_pulses_total << endl;

	// Generate pulse durations and intensities
	int N_rows = N_steps_inp;
	int N_columns = *max_element(N_pulses.begin(), N_pulses.end());

	vector<vector<double>> D, I;
	D.resize(N_rows, vector<double>(N_columns, 0.));
	I.resize(N_rows, vector<double>(N_columns, 0.));

	random_device seed_generator;
	default_random_engine generator(seed_generator());

	normal_distribution<double> D_distribution(u1, s1);
	normal_distribution<double> I_distribution(u2, s2);
	
	for (int step = 0;step < N_steps_inp;step++) {
		for (int pulse = 0;pulse < N_pulses[step];pulse++) {
			D[step][pulse] = ceil(exp(D_distribution(generator)));
			I[step][pulse] = exp(I_distribution(generator));
		}
	}
	
	// Correct for total hourly demand volume
	vector<double> corr(N_steps_inp);
	for (int step = 0;step < N_steps_inp;step++) {
		double vol_sum = 0.;
		for (int pulse = 0;pulse < N_pulses[step];pulse++) {
			vol_sum += (D[step][pulse] * I[step][pulse]);
		}
		corr[step] = (vol_sum/1000) / (flow_inp[step] * dt_h_inp * 3600);
		for (int pulse = 0;pulse < N_pulses[step];pulse++) { I[step][pulse] /= corr[step]; }
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
	for (int step = 0;step < N_steps_inp;step++) {
		for (int pulse = 0;pulse < N_pulses[step];pulse++) {

			U[step][pulse] = U_distribution(generator);
			t[step][pulse] = ((-log(U[step][pulse])) / L[step]);

			Start[k] = accumulate(t[step].begin(), t[step].end(), 0) + dt_h_inp*step * 3600;
			if (Start[k] < 0) { Start[k] == 0;}

			End[k] = Start[k] + D[step][pulse];
			if (End[k] > N_steps_inp*dt_h_inp * 3600) { End[k] = N_steps_inp*dt_h_inp * 3600; }

			Intensity[k] = I[step][pulse];
			Duration[k] = D[step][pulse];
			k++;
		}
	}

	// Sum flow from overlapping pulses
	vector<double> Q(N_steps_inp*dt_h_inp*3600,0);
	for (int pulse = 0;pulse < Start.size();pulse++) {
		for (int i = Start[pulse];i < End[pulse];++i) {
			Q[i] += Intensity[pulse];
		}
	}

	// Averaging flow rates
	int N_avg_int = dt_h_inp * 3600 / Avg_int;
	vector<double> Q_avg(N_steps_inp*N_avg_int, 0);
	int T = 0, K = 0;
	for (int step = 0;step < N_steps_inp;step++) {
		for (int interv = 0;interv < N_avg_int;interv++) {
			for (int i = T;i < T + Avg_int;i++) { Q_avg[K] += Q[i]; }
			Q_avg[K] /= (Avg_int*1000);
			T += Avg_int;
			K++;
		}
	}


	return Q_avg;
}