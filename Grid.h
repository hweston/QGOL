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
    Grid(string filename, int n);
    void GetNextState(int n, double B);
    void Collapse(int j, int k);
    vector<vector<vector<complex<double>>>> grid;

};


