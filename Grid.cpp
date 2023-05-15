//
// Created by Henry Weston on 10/11/2022.
//

#include "Grid.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <complex>

using namespace std;

Grid::Grid(string filename, int n)
        :

        grid(vector<vector<vector<complex<double>>>>
                     (n, vector<vector<complex<double>>>
                             (n,vector<complex<double>>
                                    (4,0))))


{
    time_t t;
    srand((unsigned) time(&t));
    ofstream output("outputs.txt");
    output.close();
    ifstream myfile(filename);


    if (!myfile.is_open()) {                                        /* Initialising grid with chosen seed */
        cout << "No such file" << endl;
        return;
    }
    int x = 0;
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            grid[i][j][1].real(1);                              /* Initialise all states as dead first*/
        }
    }
    while (!myfile.eof()) {
        int i, j;
        float k, l;
        myfile >> i >> j >> k >> l;


            grid[i + 1][j + 1][0].real(k);
            if(abs(pow(grid[i + 1][j + 1][0].real(),2) + pow(l,2) - 1) < 0.0001){       /* Ensures no rounding errors*/
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
            grid[i + 1][j + 1][1].real(pow((1-norm(grid[i + 1][j + 1][0])),0.5));                                   /* Inputting coefficient of 'dead' state*/




    }


    myfile.close();


}

    void Grid::GetNextState(int n, double B) {      /* Advances grid by one time step */


        ofstream myfile;
        myfile.open("outputs.txt", ios::app);
        myfile << endl;
        for (int i = 1; i < n-1; i++) {
            myfile << endl;
            for (int j = 1; j < n-1; j++) {
                if(grid[i][j][0].imag() <0){
                    myfile << "-";
                }
                myfile << abs(grid[i][j][0]) << ",";                /* Printing to output file*/

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
                    grid[i][j][3].real(pow(real(grid[i][j][0]), 2) / norm(grid[i][j][0]));      /* ratio of alive coefficient to |a| */
                    grid[i][j][3].imag(pow(imag(grid[i][j][0]), 2) / norm(grid[i][j][0]));

                    if (grid[i][j][0].imag() < 0) {
                        sign = -1;
                    }
                    else{
                    }
                }
                if (abs(neighbours[i][j]) <= 1) {                               /* Applying E operator */
                    grid[i][j][2].real(0); /* placeholder for alive state */
                    grid[i][j][1].real(1); /* dead state calculation */
                } else if ((abs(neighbours[i][j]) > 1) && (abs(neighbours[i][j]) < 2)) {
                    grid[i][j][2].real((abs(neighbours[i][j])-1)*abs(grid[i][j][0]));
                    grid[i][j][1].real(B*(2-abs(neighbours[i][j]))*(real(grid[i][j][1])+abs(grid[i][j][0])) + (abs(neighbours[i][j])-1)*real(grid[i][j][1]));
                } else if ((abs(neighbours[i][j]) >= 2) && (abs(neighbours[i][j]) <= 3)) {
                    grid[i][j][2].real(B*((3 - abs(neighbours[i][j])) * abs(grid[i][j][0])) + (abs(neighbours[i][j])-2)*(abs(grid[i][j][0]) + real(grid[i][j][1])));
                    grid[i][j][1].real(B*((3 - abs(neighbours[i][j])) * real(grid[i][j][1])));
                } else if ((abs(neighbours[i][j]) > 3) && (abs(neighbours[i][j]) <= 4)) {
                    grid[i][j][2].real(B*((4-abs(neighbours[i][j]))*(abs(grid[i][j][0]) + real(grid[i][j][1]))));
                    grid[i][j][1].real((abs(neighbours[i][j])-3)*(abs(grid[i][j][0]) + real(grid[i][j][1])));
                } else if ((abs(neighbours[i][j]) > 4)){
                    grid[i][j][2].real(0);
                    grid[i][j][1].real(1);
                }


                double norm_constant = pow(real(grid[i][j][2]),2) + pow(real(grid[i][j][1]),2);     /* Ensuring correct ratio of real/imaginary coefficients */
                grid[i][j][2].real(pow(pow(grid[i][j][2].real(),2)/norm_constant,0.5));
                grid[i][j][1].real(pow(pow(real(grid[i][j][1]),2)/norm_constant,0.5));
                if (abs(grid[i][j][0]) != 0) {
                    grid[i][j][0].real(pow(real(grid[i][j][3]) * pow(real(grid[i][j][2]),2), 0.5));
                    if(sign == -1) {
                        grid[i][j][0].imag(-pow(imag(grid[i][j][3]) * pow(real(grid[i][j][2]),2), 0.5));
                    }
                    else{
                        grid[i][j][0].imag(pow(imag(grid[i][j][3]) * pow(real(grid[i][j][2]),2), 0.5));
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
                        grid[i][j][0].real(grid[i][j][2].real()/2);

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

    void Grid::Collapse(int j, int k){          /* Measurement Process */
        int r_num = (rand() % 10000);
        grid[j][k][2].real(pow(grid[j][k][0].real(),2)/pow(grid[j][k][0].imag(),2));
        grid[j][k][3].real(pow(grid[j][k][1].real(),2)/pow(grid[j][k][1].imag(),2));

        if (10000*norm(grid[j][k][0]) > r_num){     /* Probabilistic collapse based on random number */
            if(grid[j][k][0].imag() < 0) {
                grid[j][k][0].imag(-pow((1 / (grid[j][k][2].real() + 1)), 0.5));
            }
            else{
                grid[j][k][0].imag(pow((1 / (grid[j][k][2].real() + 1)), 0.5));
            }
            grid[j][k][0].real(pow(1-pow(grid[j][k][0].imag(),2),0.5));
            grid[j][k][1] = 0;

        }
        else{
            if(grid[j][k][1].imag() < 0) {
                grid[j][k][1].imag(-pow((1 / (grid[j][k][3].real() + 1)), 0.5));
            }
            else{
                grid[j][k][1].imag(pow((1 / (grid[j][k][3].real() + 1)), 0.5));
            }
            grid[j][k][1].real(pow(1-pow(grid[j][k][1].imag(),2),0.5));
            grid[j][k][0] = 0;

        }
    }



