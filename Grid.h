//
// Created by Henry Weston on 10/11/2022.
//

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class Grid {
public:
    Grid(string filename, int n);
    void GetNextState(int n);
    vector<vector<float>> grid;

private:
    int n;

};


