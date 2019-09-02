# -*- coding: utf-8 -*-
"""
Testing WUDESIM DLL toolkit
@author: Ahmed Abokifa
"""

# Import ctypes
import ctypes

# Import numpy
import numpy

#%% Example 1 : Find dead-end branches in the network

# Import WUDESIM DLL
WUDESIM_DLL = ctypes.WinDLL('WUDESIM_LIB.dll')

# Define EPANET input and report filenames
EPANET_INP  = ctypes.c_char_p(str.encode('CTOWN.INP'))
EPANET_RPT  = ctypes.c_char_p(str.encode('CTOWN.RPT'))

# Open new EPANET project & Find deadend branches
WUDESIM_DLL.DE_ENGINE_OPEN_EPANET_PROJ(EPANET_INP,EPANET_RPT)
WUDESIM_DLL.DE_ENGINE_FIND_DEADENDS()

# Write the IDs of deadend pipes to the output file
WUDESIM_DLL.DE_WRITE_DEADEND_IDS()

# Get the number of deadend branches
N_branches = WUDESIM_DLL.DE_GET_BRAN_COUNT(0)
 
# Get the IDs of the deadend branches and pipes
bran_ID   = []
pipe_ID   = []
bran_size = []
WUDESIM_DLL.DE_GET_ID.restype = ctypes.c_char_p
for bran_idx in range(N_branches):
    bran_size.append(WUDESIM_DLL.DE_GET_BRAN_SIZE(0,bran_idx))    
    bran_ID.append(WUDESIM_DLL.DE_GET_ID(2,bran_idx,0))
    
    for pipe_idx in range(bran_size[bran_idx]):
        pipe_ID.append(WUDESIM_DLL.DE_GET_ID(0,bran_idx,pipe_idx))

#%% Example 2 : Run EPANET simulation and obtain the results

# Run EPANET simulation & Calculate the properties of the dead-end branches
WUDESIM_DLL.DE_ENGINE_RUN_EPANET_SIM()
WUDESIM_DLL.DE_ENGINE_CALC_DEADEND_PROPERTIES_EPANET()

# Write EPANET Report file and the average properties to the output file
WUDESIM_DLL.DE_WRITE_EPANET_REPORT()        
WUDESIM_DLL.DE_WRITE_DEADEND_PROPERTIES()

# Get the number of EPANET simulation steps
N_steps = WUDESIM_DLL.DE_GET_STEP_COUNT(0)

# Get the Reynolds number for all pipes as simulated by EPANET
Re_all_pipes = []
WUDESIM_DLL.DE_GET_PIPE_RESULT_EPANET.restype = ctypes.c_double
for bran_idx in range(N_branches):   

    for pipe_idx in range(bran_size[bran_idx]):
        Re_pipe_EPANET = []
        
        for step in range(N_steps):
            Re_pipe_EPANET.append(WUDESIM_DLL.DE_GET_PIPE_RESULT_EPANET(0,
                                                       bran_idx,pipe_idx,step))
    Re_all_pipes.append(Re_pipe_EPANET)


#%% Example 3 : Run WUDESIM simulation and obtain the results 

# Define WUDESIM input and report filenames
WUDESIM_INP = ctypes.c_char_p(str.encode('WUDESIM.INP'))
WUDESIM_RPT = ctypes.c_char_p(str.encode('WUDESIM.RPT'))

# Open WUDESIM project & Run WUDESIM simulation
WUDESIM_DLL.DE_ENGINE_OPEN_WUDESIM_PROJ(WUDESIM_INP,WUDESIM_RPT)
WUDESIM_DLL.DE_ENGINE_RUN_WUDESIM_SIM()

# Write WUDESIM report file
WUDESIM_DLL.DE_WRITE_WUDESIM_REPORT()

# Get WQ simulation results of EPANET and WUDESIM for a given node
bran_idx  = 0
pipe_idx  = 0
node_idx  = 0

WQ_node_EPANET  = []
WQ_node_WUDESIM = []

WUDESIM_DLL.DE_GET_NODE_RESULT_EPANET.restype  = ctypes.c_double
WUDESIM_DLL.DE_GET_NODE_RESULT_WUDESIM.restype = ctypes.c_double

for step in range(N_steps):
    WQ_node_EPANET.append(WUDESIM_DLL.DE_GET_NODE_RESULT_EPANET(0,bran_idx,
                                                                node_idx,step))

    WQ_node_WUDESIM.append(WUDESIM_DLL.DE_GET_NODE_RESULT_WUDESIM(0,bran_idx,
                                                                node_idx,step))


#%% Plot results
    
import matplotlib.pyplot as plt

# Plot water quality results
plt.rcParams.update({'font.size': 35})
fig, ax = plt.subplots(figsize=(10, 10),dpi=300)
plt.plot(range(N_steps),WQ_node_EPANET,'r-')
plt.plot(range(N_steps),WQ_node_WUDESIM,'b--')
plt.xlim([50,100])
plt.ylim([2.5,3.5])
plt.xlabel('Time (hrs)')
plt.ylabel('Chlorine Concentration (mg/L)')
plt.legend(['EPANET','WUDESIM'],frameon=False,loc='upper right')
