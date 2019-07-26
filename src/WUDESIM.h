#pragma once

#include <windows.h>

using namespace std;

#define C_DLLEXPORT extern "C" __declspec(dllexport)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Engine functions

C_DLLEXPORT int DE_RUN_FULL_SIM(const char*, const char*, const char*, const char*); //EPANET_INP, EPANET_RPT, WUDESIM_INP, WUDESIM_RPT

C_DLLEXPORT int DE_OPEN_EPANET_PROJ(const char*);  //EPANET_INP

C_DLLEXPORT int DE_FIND_DEADENDS();

C_DLLEXPORT int DE_RUN_EPANET_SIM(const char*);    //EPANET_RPT

C_DLLEXPORT int DE_CALC_DEADEND_PROPERTIES();

C_DLLEXPORT int DE_OPEN_WUDESIM_PROJ(const char*); //WUDESIM_INP

C_DLLEXPORT int DE_GENERATE_STOC_DEMAND();

C_DLLEXPORT int DE_RUN_WUDESIM_SIM(const char*);   //WUDESIM_RPT

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Define the toolkit options

enum DE_COUNT_PROPERTY {
	DE_BRAN_COUNT, // = 0
	DE_STEP_COUNT, // = 1
};

enum DE_BRAN_PROPERTY {
	DE_BRAN_SIZE, // = 0
};

enum DE_ID_PROPERTY {
	DE_PIPE_ID, // = 0
	DE_NODE_ID, // = 1
};

enum DE_VALUE_PROPERTY {
	DE_C_WUDESIM, // = 0
	DE_C_EPANET,  // = 1
	DE_Reynolds,  // = 2
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Get functions

C_DLLEXPORT int DE_GET_COUNT(const int);                                          // property idx
C_DLLEXPORT int DE_GET_BRAN_PROP(const int, const int);                           // property idx, branch idx
C_DLLEXPORT const char* DE_GET_ID(const int, const int,const int);                // property idx, branch idx, pipe idx
C_DLLEXPORT double DE_GET_PIPE_VALUE(const int, const int, const int, const int); // property idx, branch idx, pipe idx, step idx

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



