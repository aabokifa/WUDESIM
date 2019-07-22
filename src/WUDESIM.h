#pragma once

#include <windows.h>

using namespace std;

__declspec(dllexport) int DE_RUN_FULL_SIM(char*, char*, string, string); //EPANET_INP, EPANET_RPT, WUDESIM_INP, WUDESIM_RPT

__declspec(dllexport) int DE_OPEN_EPANET_PROJ(char*); //EPANET_INP

__declspec(dllexport) int DE_FIND_DEADENDS();

__declspec(dllexport) int DE_RUN_EPANET_SIM(char*); //EPANET_RPT

__declspec(dllexport) int DE_CALC_DEADEND_PROPERTIES(); 

__declspec(dllexport) int DE_OPEN_WUDESIM_PROJ(string); //WUDESIM_INP

__declspec(dllexport) int DE_GENERATE_STOC_DEMAND();

__declspec(dllexport) int DE_RUN_WUDESIM_SIM(string); //WUDESIM_RPT
