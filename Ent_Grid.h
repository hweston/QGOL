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
    Ent_Grid(string filename, string coeff_file, string coeff,int cells, vector <int> cellnum, int n, int states);
    void GetEntNextState(int n, double B, string coeff, int states, int *signal);
    /*double EntangEntropy (int states,  vector<int> cellnum, int cells);*/
    void Ent_Collapse(int n, string coeff, int *signal);
    vector<vector<vector<vector<complex<double>>>>> grid;
    vector<vector<double>> ent_states;
    vector<complex<double>> coeffs;


};
