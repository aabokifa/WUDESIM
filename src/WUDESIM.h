#pragma once

#include <windows.h>

using namespace std;

#define C_DLLEXPORT extern "C" __declspec(dllexport)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Engine functions

C_DLLEXPORT int DE_RUN_FULL_SIM(const char*, const char*, const char*, const char*); //EPANET_INP, EPANET_RPT, WUDESIM_INP, WUDESIM_RPT
C_DLLEXPORT int DE_OPEN_EPANET_PROJ(const char*);  //EPANET_INP
C_DLLEXPORT int DE_FIND_DEADENDS();
C_DLLEXPORT int DE_RUN_EPANET_SIM(const char*);    //EPANET_RPT
C_DLLEXPORT int DE_CALC_DEADEND_PROPERTIES();
C_DLLEXPORT int DE_OPEN_WUDESIM_PROJ(const char*); //WUDESIM_INP
C_DLLEXPORT int DE_GENERATE_STOC_DEMAND();
C_DLLEXPORT int DE_RUN_WUDESIM_SIM(const char*);   //WUDESIM_RPT
C_DLLEXPORT void DE_CLOSE();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Writing Functions
C_DLLEXPORT int DE_WRITE_DEADEND_IDS();
C_DLLEXPORT int DE_WRITE_DEADEND_PROPERTIES();
C_DLLEXPORT int DE_WRITE_STOCHASTIC_DEMANDS();
C_DLLEXPORT int DE_WRITE_WUDESIM_REPORT();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Get functions

C_DLLEXPORT int DE_GET_COUNT(const int);                                          // property idx
C_DLLEXPORT int DE_GET_BRAN_VALUE(const int, const int);                          // property idx, branch idx
C_DLLEXPORT const char* DE_GET_ID(const int, const int, const int);                // property idx, branch idx, pipe/node idx
C_DLLEXPORT double DE_GET_PIPE_VALUE(const int, const int, const int, const int); // property idx, branch idx, pipe idx, step idx
C_DLLEXPORT double DE_GET_NODE_VALUE(const int, const int, const int, const int); // property idx, branch idx, node idx, step idx

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Define the toolkit options

enum COUNT_OPT {
	DE_BRAN_COUNT, // = 0
	DE_EPANET_STEP_COUNT,  // = 1
	DE_WUDESIM_STEP_COUNT, // = 2
};

enum BRAN_VALUE_OPT {
	DE_BRAN_SIZE, // = 0
}; 

enum ID_OPT {
	DE_PIPE_ID, // = 0
	DE_NODE_ID, // = 1
};

enum PIPE_VALUE_OPT {
	DE_REYNOLDS,  // = 0
	DE_PECLET,    // = 1
	DE_RES_TIME,  // = 2
	DE_FLOW_EPANET,  // = 3
	DE_FLOW_WUDESIM, // = 4
	DE_LENG,      // = 5
	DE_DIAM,      // = 6
};

enum NODE_VALUE_OPT {
	DE_C_WUDESIM, // = 0
	DE_C_EPANET,  // = 1
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM log writing function

void WRITE_LOG_MSG(string);







