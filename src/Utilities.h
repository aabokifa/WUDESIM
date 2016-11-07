/*
Project:     WUDESIM ver. 1 BETA
File:        Utilities.h
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: Header file for the Utilities.cpp source file.
*/

#pragma once

#include "Classes.h"

using namespace std;

bool compare(string str1, string str2);

bool find(string str1, string str2);

vector<string> ImportFile(string);

vector<string> InputData(vector<int>&, vector<string>&, string);

void solveThomas(vector<double>, vector<double>, vector<double>, vector<double>&, int);

vector < double> interpolation(vector<double>&, vector<double>&, vector<double>&);