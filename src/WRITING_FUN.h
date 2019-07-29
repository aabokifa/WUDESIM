/*
Project:     WUDESIM ver. 1 BETA
File:        Utilities.h
Author:      Ahmed Abokifa
Date:        10/25/2016
*/

#pragma once

#include <sstream>
#include <string>

#include "Classes.h"

void write_DE_ids(Network*);
void write_DE_Properties(Network*);
void write_stoc_dems(Network*);
void write_WUDESIM_rpt(Network*);

// define template class to convert int and double to string
template<class T>
std::string toString(const T& value) {
	std::ostringstream os;
	os << value;
	return os.str();
}

