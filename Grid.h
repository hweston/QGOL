//
// Created by Henry Weston on 10/11/2022.
//

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Grid {
public:
    Grid(string filename, int n, double B, string simulation, int z);
    void GetNextState(int n, double B, string simulation,int z);
    vector<vector<vector<double>>> grid;
    string standard = "standard";
    string entangled = "entangled";

private:
    int n;

};


