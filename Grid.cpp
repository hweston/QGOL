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
#include <Eigen/Dense>
#include <complex>

using namespace std;
using namespace Eigen;

Grid::Grid(string filename, string coeff_file, string coeff,int cells, vector <int> cellnum, int n,string simulation, string etype, int z)
        :

        grid(vector<vector<vector<complex<double>>>>
                     (n, vector<vector<complex<double>>>
                             (n,vector<complex<double>>(z,0)))),
        ent_states(vector<vector<double>> (cells, vector<double>(3, 0.0))),
        coeffs(vector<complex<double>> (z))

{
    time_t t;
    srand((unsigned) time(&t));
    ofstream output("outputs.txt");
    output.close();
    ifstream myfile(filename);


    if (!myfile.is_open()) {
        cout << "No such file" << endl;
        return;
    }
    int x = 0;
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            grid[i][j][1].real(1); /* Initialise all states as dead*/
        }
    }
    while (!myfile.eof()) {
        int i, j;
        float k, l;
        myfile >> i >> j >> k >> l;

        if (simulation == "entangled") {
            for (int a = 0; a < grid[i + 1][j + 1].size() - 1; a++) {
                grid[i + 1][j + 1][a].real(k);
                if (etype == "n" && a == 1) {
                    a = grid[i + 1][j + 1].size() - 2;
                }
           }
        }
        if (simulation == "standard"){
            grid[i + 1][j + 1][0].real(k);
            if(abs(pow(grid[i + 1][j + 1][0].real(),2) + pow(l,2) - 1) < 0.0001){
                if (l<0) {
                    grid[i + 1][j + 1][0].imag(-pow(1 - pow(grid[i + 1][j + 1][0].real(), 2), 0.5));
                }
                else{
                    grid[i + 1][j + 1][0].imag(pow(1 - pow(grid[i + 1][j + 1][0].real(), 2), 0.5));
                }
            }
            else {
                grid[i + 1][j + 1][0].imag(l);
            }
            grid[i + 1][j + 1][1].real(pow((1-norm(grid[i + 1][j + 1][0])),0.5)); /* Inputting coefficient of 'dead' state*/
        }


        if (x < cells) {
            ent_states[x][0] = i + 1;
            ent_states[x][1] = j + 1;
            ent_states[x][2] = k;
            grid[ent_states[x][0]][ent_states[x][1]][z - 1].imag(0);

            x++;
        }
    }
        if (simulation == "entangled") {
            if (etype == "n") {
                for (int a = 0; a<cellnum.size(); a++){
                    for (int b = 0; b<cells; b++){
                        grid[ent_states[b][0]][ent_states[b][1]][a] = 0;
                    }
                }
                int counter = 0;
                for (int a = 0; a<cellnum.size();a++) {
                    for (int b = 0; b< cellnum[a]; b++) {
                        grid[ent_states[counter][0]][ent_states[counter][1]][a] = ent_states[counter][2];
                        counter ++;
                    }
                }

            }
            else {
                for (int b = 1; b < pow(2,z); b++) {
                    for (int a = 0; a < z; a++) {
                        grid[ent_states[a][0]][ent_states[a][1]][b - 1] = (b & (1 << a)) >> a;
                    }
                }
            }
        }

    myfile.close();

    if (coeff == "y"){
        ifstream cfile(coeff_file);
        if (!cfile.is_open()) {
            cout << "No such file" << endl;
            return;
        }
        int k = 0;
        while (!cfile.eof()) {
            double i, j;
            cfile >> i >> j;

            complex<double> ctemp (i,j);
            coeffs[k] = ctemp;
            /*norm += abs(ctemp); */
            k++;

        }
        cfile.close();
    }
}

    void Grid::GetNextState(int n, double B, string simulation,string coeff, int z, string etype, int *signal) {

        if (simulation == "standard") {
            ofstream myfile;
            myfile.open("outputs.txt", ios::app);
            myfile << endl;
            for (int i = 1; i < n-1; i++) {
                myfile << endl;
                for (int j = 1; j < n-1; j++) {
                    myfile << abs(grid[i][j][0]) << ',';

                }

            }
            myfile.close();

            vector<vector<complex<double>>> neighbours(n, vector<complex<double>>(n, 0.0));


            for (int i = 1; i < n-1; i++) {                   // Checking number of neighbours
                for (int j = 1; j < n-1; j++) {


                    for (int a = -1; a < 2; a++) {              //Counting number of neighbours for each cell
                        for (int b = -1; b < 2; b++) {
                            if ((a == 0) && (b == 0)) {
                            } else {
                                neighbours[i][j] += grid[i + a][j + b][0];
                            }
                        }
                    }
                }
            }
            for (int i = 1; i < n-1; i++) {                   // Checking number of neighbours and assigning new state
                for (int j = 1; j < n-1; j++) {
                    int sign = 0;
                    if (abs(grid[i][j][0]) != 0) {
                        grid[i][j][3].real(pow(real(grid[i][j][0]), 2) / norm(grid[i][j][0]));
                        grid[i][j][3].imag(pow(imag(grid[i][j][0]), 2) / norm(grid[i][j][0]));

                        if (grid[i][j][0].imag() < 0) {
                            sign = -1;
                        }
                        else{
                        }
                    }
                    if (abs(neighbours[i][j]) <= 1) {
                        grid[i][j][2].real(0); /* placeholder for alive state */
                        grid[i][j][1].real(1); /* dead state calculation */
                    } else if ((abs(neighbours[i][j]) > 1) && (abs(neighbours[i][j]) <= 2)) {
                        grid[i][j][2].real((abs(neighbours[i][j])-1)*abs(grid[i][j][0]));
                        grid[i][j][1].real(B*(2-abs(neighbours[i][j]))*(real(grid[i][j][1])+abs(grid[i][j][0])) + (abs(neighbours[i][j])-1)*real(grid[i][j][1]));
                    } else if ((abs(neighbours[i][j]) > 2) && (abs(neighbours[i][j]) <= 3)) {
                        grid[i][j][2].real(B*((3 - abs(neighbours[i][j])) * abs(grid[i][j][0])) + (abs(neighbours[i][j])-2)*(abs(grid[i][j][0]) + real(grid[i][j][1])));
                        grid[i][j][1].real(B*((3 - abs(neighbours[i][j])) * real(grid[i][j][1])));
                    } else if ((abs(neighbours[i][j]) > 3) && (abs(neighbours[i][j]) <= 4)) {
                        grid[i][j][2].real(B*((4-abs(neighbours[i][j]))*(abs(grid[i][j][0]) + real(grid[i][j][1]))));
                        grid[i][j][1].real((abs(neighbours[i][j])-3)*(abs(grid[i][j][0]) + real(grid[i][j][1])));
                    } else if ((abs(neighbours[i][j]) > 4)){
                        grid[i][j][2].real(0);
                        grid[i][j][1].real(1);
                    }

                    double norm_constant = real(grid[i][j][2]) + real(grid[i][j][1]);
                    grid[i][j][2].real(grid[i][j][2].real()/norm_constant);
                    grid[i][j][1].real(real(grid[i][j][1])/norm_constant);
                    if (abs(grid[i][j][0]) != 0) {
                        grid[i][j][0].real(pow(real(grid[i][j][3]) * real(grid[i][j][2]), 0.5));
                        if(sign == -1) {
                            grid[i][j][0].imag(-pow(imag(grid[i][j][3]) * real(grid[i][j][2]), 0.5));
                        }
                        else{
                            grid[i][j][0].imag(pow(imag(grid[i][j][3]) * real(grid[i][j][2]), 0.5));
                        }
                    }
                    else{
                        if (neighbours[i][j].real() == 0){
                            if (neighbours[i][j].imag()<0) {
                                grid[i][j][0].imag(-grid[i][j][2].real());
                            }
                            else{
                                grid[i][j][0].imag(grid[i][j][2].real());
                            }
                        }
                        else if(neighbours[i][j].imag() == 0){
                            grid[i][j][0].real(grid[i][j][2].real());
                        }
                        else{
                            grid[i][j][0].real(grid[i][j][2].real()*pow(2,0.5)/2);

                            if (neighbours[i][j].imag()<0) {
                                grid[i][j][0].imag(-grid[i][j][0].real());
                            }
                            else{
                                grid[i][j][0].imag(grid[i][j][0].real());
                            }
                        }
                    }

                    grid[i][j][2].real(0);
                }
            }


        }



        /*
        Getting next state for entangled states
        */
        int states = pow(2, z);
        if (etype == "n") {
            states = 2;
        }
        int prob = 0;
        if (simulation == "entangled") {
            vector<vector<double>> wavefunction(n, vector<double>(n, 0.0));
            ofstream myfile;
            myfile.open("outputs.txt", ios::app);
            myfile << endl;

            for (int i = 0; i < n - 2; i++) {
                myfile << endl;
                for (int k = 0; k < states; k++) {
                    for (int j = 0; j < n - 2; j++) {

                        myfile << grid[i + 1][j + 1][k] << ",";
                        wavefunction[i + 1][j + 1] += grid[i + 1][j + 1][k];

                    }
                    myfile << "\t";

                }


                for (int j = 0; j < n - 2; j++) {
                    if (*signal == 2) {
                        myfile << wavefunction[i + 1][j + 1] << ",";
                    } else {
                        myfile << wavefunction[i + 1][j + 1] / states << ",";
                    }
                }
            }


            myfile.close();

            vector<vector<vector<double>>> neighbours(n, vector<vector<double>>(n, vector<double>(pow(2, z), 0.0)));

            for (int k = 0; k < (pow(2, z)); k++) {
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
                for (int i = 1;i < n - 1; i++) {                   // Checking number of neighbours and assigning new state
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
                        if (j == 0 || j == n-2) {
                            if (grid[i][j][k] != 0 && *signal == 0) {
                                if (coeff == "y"){
                                    vector<int> ranges(2);
                                    double ctotal;
                                    for (int m = 0; m<2;m++){
                                        ctotal += (abs(coeffs[m])/norm);
                                        ranges[m] = 10000*ctotal;
                                    }
                                    int r_num = (rand() % 10000);
                                    int m = 0;
                                    while(r_num>ranges[m]){
                                        m++;
                                        prob++;
                                    }

                                }
                                else {
                                    prob = (rand() % states);
                                }
                                *signal = 1;


                            }
                        }
                    }
                }


            }
            if (*signal == 1) {
                for (int l = 0; l < n - 2; l++) {
                    for (int m = 0; m < states; m++) {
                        for (int p = 0; p < n - 2; p++) {
                            if (m != prob) {
                                grid[l + 1][p + 1][m] = 0;

                            }
                        }
                    }
                }
                *signal = 2;
            }

        }
    }
/*
    double Grid::EntangEntropy(int states, vector<int> cellnum, int cells) {

        int maxcells = *max_element(cellnum.begin(), cellnum.end());
        vector<vector<int>> cellstates(states, vector<int>(maxcells));
        int counter = 0;
        for (int j = 0; j<states; j++){
          for (int i = 0; i < cellnum[j]; i++) {
                cellstates[j][i] = grid[ent_states[counter][0]][ent_states[counter][1]][j];
                counter ++;
            }
        }

        MatrixXd density_matrix(num,num);
        density_matrix = MatrixXd::Zero(num,num);
        for (int j = 0; j<states;j++){
            VectorXd phi1 = VectorXd::Map(cells[j].data(), cells[j].size());
            for (int i = 0; i<states; i++) {
                VectorXd phi2 = VectorXd::Map(cells[i].data(), cells[i].size());
                density_matrix += phi1 * phi2.transpose();
            }
        }

        cout << density_matrix;

    }
*/
        void Grid::Collapse(int j, int k){
        int r_num = (rand() % 10000);
        if (10000*norm(grid[j][k][0]) > r_num){
            grid[j][k][0] = 1;
            grid[j][k][1] = 0;
        }
        else{
            grid[j][k][0] = 0;
            grid[j][k][1] = 1;
        }

}

