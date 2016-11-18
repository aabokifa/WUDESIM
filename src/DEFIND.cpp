/*
Project:     WUDESIM ver. 1 BETA
File:        DEFIND.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: This function studies the connectivity of each junction in the network to determine if it is a dead-end, 
			 finds all dead-end branches in the network, and generates a list of dead-end pipe IDs to the file “DEpipeID.txt”.
*/

#include <iostream> 
#include <fstream>  
#include <vector>  
#include <map>
#include <algorithm>
#include <sstream>  
#include <string>  
#include <iterator>
#include <numeric>
#include <stdio.h>
#include <tchar.h>
#include <SDKDDKVer.h>

#include "Classes.h"
#include "WUDESIMmain.h"
#include "Utilities.h"

using namespace std;


int DEFIND(Network* net)
{
	//Find number of connections for node_1 and node_2 of each pipe
	int N_pipes = net->pipes.size();

	vector<string> node_1_2(2 * N_pipes);
	int j = 0;
	for (int i = 0;i < N_pipes;++i) {
		node_1_2[j] = net->pipes[i].node_1;
		j++;
		node_1_2[j] = net->pipes[i].node_2;
		j++;
	}

	for (int i = 0;i < N_pipes;i++) {

		//Find number of connections for node_1 and node_2 of each pipe
		for (int j = 0;j < node_1_2.size();++j) {
			if (compare(net->pipes[i].node_1, node_1_2[j])) {
				net->pipes[i].node_1_conn++;
			}
			if (compare(net->pipes[i].node_2, node_1_2[j])) {
				net->pipes[i].node_2_conn++;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Find dead end branches

	int branch = 0;        // First DE branch
	int N_branches = 0;

	string dead_node, other_node;

	int other_node_conn;

	for (int i = 0;i < N_pipes;i++) {

		// Condition 1: Check if one of the two nodes is connected to only one pipe
		if (net->pipes[i].node_1_conn == 1) {

			dead_node = net->pipes[i].node_1;
			other_node = net->pipes[i].node_2;
			other_node_conn = net->pipes[i].node_2_conn; //Number of connections to the other node
		}
		else if (net->pipes[i].node_2_conn == 1) {

			dead_node = net->pipes[i].node_2;
			other_node = net->pipes[i].node_1;
			other_node_conn = net->pipes[i].node_1_conn; //Number of connections to the other node
		}

		else { goto nextBranch; } //if it's not a dead-end then go directly to the next pipe


		//Condition 2: Check DE node is not a source
		int N_sources = net->demand_sources.size();
		for (int j = 0;j < N_sources;++j) {
			if (compare(dead_node, net->demand_sources[j])) { goto nextBranch; }
		}

		// Condition 3: Check DE node is not a tank
		int N_tanks = net->tanks.size();
		for (int j = 0;j < N_tanks;++j) {
			if (compare(dead_node, net->tanks[j])) { goto nextBranch; }
		}

		// Condition 4: Check DE node is not a reservoir
		int N_reserv = net->reservoirs.size();
		for (int j = 0;j < N_reserv;++j) {
			if (compare(dead_node, net->reservoirs[j])) { goto nextBranch; }
		}

		// Condition 5: Check DE node is not a pump
		int N_pumps = net->pumps.size();
		for (int j = 0;j < N_pumps;++j) {
			if (compare(dead_node, net->pumps[j].start)) { goto nextBranch; }
			if (compare(dead_node, net->pumps[j].end)) { goto nextBranch; }
		}

		// Condition 6: Check DE node is not a valve
		int N_valves = net->valves.size();

		for (int j = 0;j < N_valves;++j) {
			if (compare(dead_node, net->valves[j].start)) { goto nextBranch; }
			if (compare(dead_node, net->valves[j].end)) { goto nextBranch; }
		}

		// if you come here then you are a dead-end!
		N_branches++;
		net->DE_branches.resize(N_branches);
		net->DE_branches[branch].pipe_id.push_back(net->pipes[i].id);
		net->DE_branches[branch].pipe_index.push_back(i);
		net->DE_branches[branch].branch_size++;

		//Now check to see if preceding pipes are to be added to the branch
		if (other_node_conn == 2) {
			
			int k = i;

		anotherPipe:; //Will come back to here if more pipes are to be added

			// Check that the upstream node is not a quality source
			for (int source = 0; source < net->quality_sources.size();source++) {
				if (compare(other_node, net->quality_sources[source])) {branch++;goto nextBranch;}
			}

			for (int p = 0;p < N_pipes;p++) {

				if (compare(net->pipes[p].node_1, other_node) && p != k) {

					net->DE_branches[branch].pipe_id.push_back(net->pipes[p].id);
					net->DE_branches[branch].pipe_index.push_back(p);
					net->DE_branches[branch].branch_size++;

					if (net->pipes[p].node_2_conn == 2) {
						other_node = net->pipes[p].node_2;
						k = p;
						goto anotherPipe;
					}
					else { branch++; goto nextBranch; }
				}

				else if (compare(net->pipes[p].node_2, other_node) && p != k) {

					net->DE_branches[branch].pipe_id.push_back(net->pipes[p].id);
					net->DE_branches[branch].pipe_index.push_back(p);
					net->DE_branches[branch].branch_size++;


					if (net->pipes[p].node_1_conn == 2) {
						other_node = net->pipes[p].node_1;
						k = p;
						goto anotherPipe;
					}
					else { branch++; goto nextBranch; }


				}
			}

		}
		else { branch++; goto nextBranch; }

	nextBranch:;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Write dead end branch ids to DEpipeID.txt file
	ofstream ofs;
	ofs.open("DEpipeID.txt", ios::out | ios::trunc);
	if (ofs.is_open()) {
		ofs << "Deadend pipe IDs:" << endl;

		for (int branch = 0;branch < N_branches;branch++) {
			ofs << "DE Branch No. " << branch + 1 << " :" << endl;

			for (int pipe = 0;pipe < net->DE_branches[branch].branch_size;pipe++) {
				ofs << net->DE_branches[branch].pipe_id[pipe] << endl;

			}
		}
		ofs.close();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	


	return N_branches;
}