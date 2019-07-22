/*
Project:     WUDESIM ver. 1 BETA
File:        WUDESIMmain.h
Author:      Ahmed Abokifa
Date:        10/25/2016
*/


#pragma once


#include "Classes.h"
#include <vector> 

using namespace std;

int OP_EPANET_INP(Network*);

int FIND_DE_BRANCHES(Network*);

int OP_WUDESIM_INP(Network*);

int RUN_EPANET_SIM(Network*);

int CALC_CORR_FACT(Network*);

int CALC_DE_PROPERTIES(Network*);

int GEN_STOC_DEM(Network*);

int RUN_WUDESIM_SIM(Network*);
