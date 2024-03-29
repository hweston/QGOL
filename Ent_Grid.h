//
// Created by Henry Weston on 10/11/2022.
//

#include <iostream>
#include <vector>
#include <string>
#include <complex>

using namespace std;

class Ent_Grid {
public:
    Ent_Grid(string filename, string coeff_file,int cells, vector <int> cellnum, int n, int states);
    void GetEntNextState(int n, double B, int states, int *signal);
    double EntangleEntropy (int n, int states, string ent_file, int sub_cells);
    void Ent_Collapse(int n, int *signal, string coeff_file);
    vector<vector<vector<vector<complex<double>>>>> grid;
    vector<vector<double>> ent_states;
    vector<complex<double>> coeffs;


};
