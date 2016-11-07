/*
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.h
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: Header file for WUDESIMmain.cpp source file.
*/


#pragma once


#include "Classes.h"


using namespace std;

int OpenEPANETinp(string, Network*);

int DEFIND(Network*);

int XJFIND(Network*);

int OpenEPANETrpt(string, Network*);

int WQSIM(string, Network*);

vector<double> DEMGEN(vector<double>& flow_inp, double dt_q, double dt_h, vector<double>&);
