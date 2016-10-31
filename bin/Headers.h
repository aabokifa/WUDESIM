/*
File:   Headers.h
Author: Ahmed Abokifa
Date:   10/25/2016
Desc:	This header file contains the declarations for all functions imbedded in the different source codes of the program.
*/

#pragma once
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

extern ofstream logfile;
extern ofstream myrptfile;

bool compare(string str1, string str2);

bool find(string str1, string str2);

vector<string> ImportFile(string);

vector<string> InputData(vector<int>&, vector<string>&, string);

void solveThomas(vector<double>, vector<double>, vector<double>, vector<double>&, int); 

vector < double> interpolation(vector<double>&, vector<double>&, vector<double>&);

int WUDESIM(int, vector<string>&);

vector<double> DEMGEN(vector<double>& flow_inp, double dt_q, double dt_h, vector<double>&);

