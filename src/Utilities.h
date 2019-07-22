/*
Project:     WUDESIM ver. 1 BETA
File:        Utilities.h
Author:      Ahmed Abokifa
Date:        10/25/2016
*/

#pragma once

#include "Classes.h"
#include <vector>

using namespace std;

bool compare_str(string str1, string str2);

bool find_str(string str1, string str2);

vector<string> ImportFile(string);

vector<string> InputData(vector<int>&, vector<string>&, string);

void solveThomas(vector<double>, vector<double>, vector<double>, vector<double>&, int);

vector < double> interpolation(vector<double>&, vector<double>&, vector<double>&);

double avrg(vector<double>&);