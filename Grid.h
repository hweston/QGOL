//
// Created by Henry Weston on 10/11/2022.
//

#include <iostream>
#include <vector>
#include <string>
#include <complex>

using namespace std;

class Grid {
public:
    Grid(string filename, string coeff_file, string coeff,int cells, vector <int> cellnum, int n, string simulation, string etype, int z);
    void GetNextState(int n, double B, string simulation, string coeff, int z, string etype, int* signal);
    void Collapse(int j, int k);
    double EntangEntropy (int states,  vector<int> cellnum, int cells);
    vector<vector<vector<complex<double>>>> grid;
    vector<vector<double>> ent_states;
    vector<complex<double>> coeffs;
    double norm;

};


