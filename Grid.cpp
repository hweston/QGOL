//
// Created by Henry Weston on 10/11/2022.
//

#include "Grid.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>

using namespace std;

Grid::Grid(string filename, int n, double B,string simulation, int z)
:
    n(n),
    grid(vector<vector<vector<double>>>
            (n, vector<vector<double>>
                                            (n,vector<double>(pow(2,z), 0.0))))
{

    ofstream output("outputs.txt");
    output.close();
    ifstream myfile (filename);


    if (!myfile.is_open()) {
        cout << "No such file" << endl;
        return;
    }
    int x = 0;
    vector<vector<double>> ent_states(n, vector<double>(n, 0.0));
    while(!myfile.eof()) {
        int i, j;
        float k;
        myfile >> i >> j >> k;

        for (int a = 0; a < grid[i+1][j+1].size(); a++) {
            grid[i + 1][j + 1][a] = k;
        }


        if (x < z) {
            ent_states[x][0] = i + 1;
            ent_states[x][1] = j + 1;
            grid[ent_states[x][0]][ent_states[x][1]][pow(2, z) - 1] = 0;
            x++;
        }


    }

    if (simulation == entangled){
        for(int b=1;b<pow(2,z);b++){
            for(int a = 0;a<z;a++) {
                grid[ent_states[a][0]][ent_states[a][1]][b - 1] = (b & (1 << a)) >> a;
                    }

                }

        }
    myfile.close();
}

void Grid::GetNextState(int n, double B, string simulation, int z) {

    if (simulation == standard) {
        ofstream myfile;
        myfile.open("outputs.txt");
        myfile << endl;

        for(int i =1;i<n-1; i++){
            myfile << endl;
            for(int j = 1;j<n-1;j++){
                myfile << grid[i][j][0] << ',';

            }

        }
        myfile.close();

            vector<vector<double>> neighbours(n, vector<double>(n, 0.0));
            vector<vector<double>> growth_rate(n, vector<double>(n, 0.0));

        for (int i = 1; i < n - 1; i++) {                   // Checking number of neighbours
                for (int j = 1; j < n - 1; j++) {


                    for (int a = -1; a < 2; a++) {              //Counting number of neighbours for each cell
                        for (int b = -1; b < 2; b++) {
                            if ((a == 0) && (b == 0)) {
                            } else {
                                neighbours[i][j] += (double) grid[i + a][j + b][0];
                            }
                        }
                    }
                }
            }
            for (int i = 1; i < n - 1; i++) {                   // Checking number of neighbours and assigning new state
                for (int j = 1; j < n - 1; j++) {
                    if (neighbours[i][j] < 2) {
                        growth_rate[i][j] = (double) exp(B * (neighbours[i][j] - 2));
                        grid[i][j][0] = (double) grid[i][j][0] *
                                     (growth_rate[i][j] - (growth_rate[i][j] - 0.9) * sin(M_PI * grid[i][j][0] / 2));
                    } else if ((neighbours[i][j] >= 2) && (neighbours[i][j] <= 3)) {
                        growth_rate[i][j] = (double) (3 / 2) * (sin(M_PI * (neighbours[i][j] - 2)) + 1);
                        grid[i][j][0] = (double) 0.1 + grid[i][j][0] * (growth_rate[i][j] - (growth_rate[i][j] - 0.9) *
                                                                                      sin(M_PI * grid[i][j][0] / 2));
                    } else if (neighbours[i][j] > 3) {
                        growth_rate[i][j] = (double) exp(-B * (neighbours[i][j] - 3));
                        grid[i][j][0] = (double) grid[i][j][0] *
                                     (growth_rate[i][j] - (growth_rate[i][j] - 0.9) * sin(M_PI * grid[i][j][0] / 2));
                    }


                }
            }
        }
        if (simulation == entangled) {
            vector<vector<double>> wavefunction(n, vector<double>(n, 0.0));
            ofstream myfile;
            myfile.open("outputs.txt", ios::app);
            myfile << endl;

            for (int i = 0; i < n - 2; i++) {
                myfile << endl;
                for (int k = 0; k < pow(2,z); k++) {
                    for (int j = 0; j < n - 2; j++) {
                        myfile << grid[i+1][j+1][k] << ',';
                        wavefunction[i+1][j+1]+=grid[i+1][j+1][k];

                    }
                    myfile << "\t";
                }

                for (int j = 0; j < n - 2; j++) {
                    myfile << wavefunction[i+1][j+1]/pow(2,z) << ",";
                }
            }






            myfile.close();
            vector<vector<vector<double>>>neighbours(n, vector<vector<double>>(n,vector<double>(pow(2,n), 0.0)));

            for (int k = 0; k<(pow(2,z)); k++) {
                for (int i = 1; i < n - 1; i++) {                   // Checking number of neighbours
                    for (int j = 1; j < n - 1; j++) {


                        for (int a = -1; a < 2; a++) {              //Counting number of neighbours for each cell
                            for (int b = -1; b < 2; b++) {
                                if ((a == 0) && (b == 0)) {
                                } else {
                                    neighbours[i][j][k] += grid[i + a][j + b][k];
                                }
                            }
                        }
                    }
                }
                for (int i = 1;
                     i < n - 1; i++) {                   // Checking number of neighbours and assigning new state
                    for (int j = 1; j < n - 1; j++) {
                        if ((neighbours[i][j][k] < 2) || (neighbours[i][j][k] > 3)) {
                            grid[i][j][k] = 0;

                        } else if (neighbours[i][j][k] == 2) {
                            if (grid[i][j][k] == 1) {
                                grid[i][j][k] = 1;
                            } else if (grid[i][j][k] == 0) {
                                grid[i][j][k] = 0;
                            }
                        } else if (neighbours[i][j][k] == 3) {
                            grid[i][j][k] = 1;
                        }

                    }
                }
            }
        }
        }




