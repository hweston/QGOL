//
// Created by Henry Weston on 10/11/2022.
//

#include "Grid.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

Grid::Grid(string filename, int n)
:
    n(n),
    grid(vector<vector<float>>(n, vector<float>(n,0)))
{
    time_t t;
    srand((unsigned) time(&t));
    ofstream output("outputs.txt");
    output.close();
    ifstream myfile;
    myfile.open(filename);

    if (!myfile.is_open()) {
        cout << "No such file" << endl;
        return;
    }

    while(!myfile.eof()){
        int i, j;
        float k;
        myfile >> i >> j >> k;
        grid[i+1][j+1] = k;
    }

    myfile.close();

}

void Grid::GetNextState(int n) {
    {
        ofstream myfile;
        myfile.open("outputs.txt", ios::app);
        myfile << endl;

        for(int i =1;i<n-1; i++){
            myfile << endl;
            for(int j = 1;j<n-1;j++){
                myfile << grid[i][j];
            }

        }
        myfile.close();

        vector<vector<int>> neighbours(n, vector<int>(n, 0));
        for (int i = 1; i < n - 1; i++) {                   // Checking number of neighbours
            for (int j = 1; j < n - 1; j++) {


                for (int a = -1; a < 2; a++) {              //Counting number of neighbours for each cell
                    for (int b = -1; b < 2; b++) {
                        if ((a == 0) && (b == 0)) {
                        } else {
                            neighbours[i][j] += grid[i + a][j + b];
                        }
                    }
                }
            }
        }
        for (int i = 1; i < n - 1; i++) {                   // Checking number of neighbours and assigning new state
            for (int j = 1; j < n - 1; j++) {
                if ((neighbours[i][j] < 2) || (neighbours[i][j] > 3)) {
                    int prob = rand() % 1000;
                    if (prob <= 997) {
                        grid[i][j] = 0;
                    }
                    else
                        grid[i][j] = 1;
                } else if (neighbours[i][j] == 2) {
                    if (grid[i][j] == 1) {
                        grid[i][j] = 1;
                    } else if (grid[i][j] == 0) {
                        grid[i][j] = 0;
                    }
                } else if (neighbours[i][j] == 3) {
                    grid[i][j] = 1;
                }

            }
        }

    }



}