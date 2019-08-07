#pragma once

#include <windows.h>

using namespace std;

#define C_DLLEXPORT extern "C" __declspec(dllexport)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Run full simulation function

C_DLLEXPORT int DE_RUN_FULL_SIM(const char*, const char*, const char*, const char*); //EPANET_INP, EPANET_RPT, WUDESIM_INP, WUDESIM_RPT

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Engine functions

C_DLLEXPORT int DE_ENGINE_OPEN_EPANET_PROJ(const char*, const char*);  //EPANET_INP, EPANET_RPT
C_DLLEXPORT int DE_ENGINE_FIND_DEADENDS();
C_DLLEXPORT int DE_ENGINE_RUN_EPANET_SIM();    
C_DLLEXPORT int DE_ENGINE_CALC_DEADEND_PROPERTIES_EPANET();
C_DLLEXPORT int DE_ENGINE_OPEN_WUDESIM_PROJ(const char*, const char*); //WUDESIM_INP, WUDESIM_RPT
C_DLLEXPORT int DE_ENGINE_GENERATE_STOC_DEMAND();
C_DLLEXPORT int DE_ENGINE_RUN_WUDESIM_SIM();  
C_DLLEXPORT int DE_ENGINE_CLOSE_WUDESIM_PROJ();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Writing Functions

C_DLLEXPORT int DE_WRITE_DEADEND_IDS();
C_DLLEXPORT int DE_WRITE_DEADEND_PROPERTIES();
C_DLLEXPORT int DE_WRITE_STOCHASTIC_DEMANDS();
C_DLLEXPORT int DE_WRITE_WUDESIM_REPORT();
C_DLLEXPORT int DE_WRITE_EPANET_REPORT();


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM Get functions

C_DLLEXPORT int DE_GET_COUNT(const int);                                           // property idx
C_DLLEXPORT int DE_GET_BRAN_PROPERTY(const int, const int);                        // property idx, branch idx
C_DLLEXPORT double DE_GET_PIPE_PROPERTY(const int, const int, const int);          // property idx, branch idx, pipe idx
C_DLLEXPORT double DE_GET_PIPE_RESULT_EPANET(const int, const int, const int, const int); // property idx, branch idx, pipe idx, step idx
C_DLLEXPORT double DE_GET_PIPE_RESULT_WUDESIM(const int, const int, const int, const int); // property idx, branch idx, pipe idx, step idx
C_DLLEXPORT double DE_GET_NODE_RESULT_EPANET(const int, const int, const int, const int); // property idx, branch idx, node idx, step idx
C_DLLEXPORT double DE_GET_NODE_RESULT_WUDESIM(const int, const int, const int, const int); // property idx, branch idx, node idx, step idx
C_DLLEXPORT double DE_GET_STOC_FLOW(const int, const int, const int, const int);   // property idx, branch idx, pipe idx, step idx
C_DLLEXPORT const char* DE_GET_ID(const int, const int, const int);                // property idx, branch idx, pipe/node idx

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Define the toolkit options

enum COUNT_OPT {
	DE_BRAN_COUNT,            // = 0
	DE_EPANET_STEP_COUNT,     // = 1
	DE_STOCHASTIC_STEP_COUNT, // = 2
};

enum BRAN_VALUE_OPT {
	DE_BRAN_SIZE, // = 0
}; 

enum ID_OPT {
	DE_PIPE_ID, // = 0
	DE_NODE_ID, // = 1
	DE_BRAN_ID, // = 2
};

enum PIPE_PROPERTY_OPT {
	DE_LENGTH,           // = 0
	DE_DIAMETER,         // = 1
};

enum PIPE_RESULT_EPANET_OPT {
	DE_REYNOLDS_EPANET,     // = 0
	DE_RES_TIME_EPANET,     // = 1
	DE_FLOW_EPANET,         // = 2
};

enum NODE_RESULT_EPANET_OPT {
	DE_QUAL_EPANET,      // = 0
	DE_DEMAND_EPANET, // = 1
};

enum PIPE_RESULT_WUDESIM_OPT {
	DE_REYNOLDS_WUDESIM,     // = 0
	DE_RES_TIME_WUDESIM,     // = 1
	DE_PECLET_WUDESIM,       // = 2
};

enum NODE_RESULT_WUDESIM_OPT {
	DE_QUAL_WUDESIM,       // = 0
};

enum STOC_FLOW_OPT {
	DE_FLOW_STOCHASTIC,   // = 0
	DE_DEMAND_STOCHASTIC, // = 1
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WUDESIM log writing function

void WRITE_LOG_MSG(string);







