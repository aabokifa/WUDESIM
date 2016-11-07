/*
Project:     WUDESIM ver. 1 BETA
File:        Utilities.cpp
Author:      Ahmed Abokifa
Date:        10/25/2016
Description: This source file contains multiple small functions that are frequently utilized during the execution of other source codes.
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

bool compare(string str1, string str2) //Compare String 1 to the begining of String 2
{
	if (str2.compare(0, str1.length() + 10, str1) == 0) return true;
	else return false;
}

bool find(string str1, string str2) //Find String 1 in String 2 case insensitive
{
	char c1;
	char c2=*str1.begin();   //c2 is the first character in str1
	char c3;
	char c4;

	for (string::iterator it1 = str2.begin();it1 != str2.end()-1;++it1) { //Read str2 characters one by one
		
		c1 = *it1; //c1 is the current search character in str2
		
		if (towlower(c1)==towlower(c2)) { //Compare c2 with c1 (the first character in str1)
			
			string::iterator it2 = it1; //mark the iterator location where a match happened

			for (string::iterator it3 = str1.begin();it3 != str1.end();++it3) { //Compare the rest of str2 characters with the rest of str1 characters
				c3 = *it2; //c3 is the search character in str2
				c4 = *it3; //c4 is the search character in str1
				if (towlower(c3) == towlower(c4)) { ++it2;  } //if c3 and c4 matched, move foreward one character in both strings
				else { goto nextmatch; }
			}
			return true; //if all characters in str1 matched with str2 then true
		}
	nextmatch:;
	}
	return false; //no match was found
}

vector<string> ImportFile(string FileName)
{
	string line;
	vector<string> Import;
	ifstream InpFile(FileName);
	
	if (InpFile.is_open()) {
		while (getline(InpFile, line))
		{
			if (line.length() != 0 && line[0] != ';') { Import.push_back(line); }
		}
		return Import;
	}
	else
	{
		return{};
	}
}

vector<string> InputData(vector<int>& index, vector<string>& A1, string Header)
{
	vector<int> k_start;
	vector<int> k_end;
	for (int i = 0;i < index.size();++i) {
		if (find(Header,A1[index[i]])) {
			k_start.push_back(index[i] + 1);
			k_end.push_back(index[i + 1]);
		}
	}
	vector<string> data;
	for (int i = 0;i < k_start.size();++i) {
		data.insert(data.end(), A1.begin() + k_start[i], A1.begin() + k_end[i]);
	}
	return data;
}

void solveThomas(vector<double> a, vector<double> b, vector<double> c, vector<double> &d, int n) { //Solve Tridiagonal matrix with thomas algorithm
	n--; // since we start from x0 (not x1) 
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

	for (int i = n; i-- > 0;) {
		d[i] -= c[i] * d[i + 1];
	}
}

vector < double> interpolation(vector<double>& X, vector<double>& Y, vector<double>& X_new) {
	int N = X.size();
	int N_new = X_new.size();
	vector<double> Y_new(N_new);
	for (int j = 0;j < N_new;++j) {
		for (int i = 0;i < N;++i) {
			if (abs(X_new[j] - X[i]) / X_new[j] <= 1E-6) { X_new[j] = X[i]; break; }
			else if ((X_new[j] >= X[i]) && (X[i + 1] >= X_new[j])) { Y_new[j] = Y[i] + (Y[i + 1] - Y[i])*(X_new[j] - X[i]) / (X[i + 1] - X[i]); break; }
		}
	}
	return Y_new;
}