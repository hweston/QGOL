//
// Created by Henry Weston on 10/11/2022.
//

#include "Ent_Grid.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include <complex>
#include <unsupported/Eigen/MatrixFunctions>


using namespace std;
using namespace Eigen;

Ent_Grid::Ent_Grid(string filename, string coeff_file,int cells, vector <int> cellnum, int n, int states)
        :

        grid(vector<vector<vector<vector<complex<double>>>>>
                     (n, vector<vector<vector<complex<double>>>>
                             (n,vector<vector<complex<double>>>
                                     (states+1, vector<complex<double>>
                                             (4,0))))),
        ent_states(vector<vector<double>> (cells, vector<double>(3, 0.0))),
        coeffs(vector<complex<double>> (states))

{
    time_t t;
    srand((unsigned) time(&t));
    ofstream output("outputs.txt");
    output.close();
    ifstream myfile(filename);


    if (!myfile.is_open()) {                            /* Initialising grid with chosen seed */
        cout << "No such file" << endl;
        return;
    }
    int x = 0;
    for (int i = 0; i<n; i++){
        for (int j = 0; j<n; j++){
            for (int k = 0; k<states; k++)
            grid[i][j][k][1].real(1);       /* Initialise all states as dead*/
        }
    }
    while (!myfile.eof()) {
        int i, j;
        float k, l;
        myfile >> i >> j >> k >> l;


        if (x < cells) {                        /* Assigning cells their state */
            ent_states[x][0] = i + 1;
            ent_states[x][1] = j + 1;
            ent_states[x][2] = k;
            grid[ent_states[x][0]][ent_states[x][1]][states - 1][1].real(0);

            x++;
        }
    }


    for (int a = 0; a<states; a++){
        for (int b = 0; b<cells; b++){
            grid[ent_states[b][0]][ent_states[b][1]][a][1].real(1);         /* Initialising entangled cells as dead - don't want overlap between states */
        }
    }
    int counter = 0;
    for (int a = 0; a<states;a++) {
        for (int b = 0; b< cellnum[a]; b++) {                                    /* Re-initialising chosen entangled cells as dead */
            grid[ent_states[counter][0]][ent_states[counter][1]][a][0].real(ent_states[counter][2]);
            counter ++;
        }
    }

    myfile.close();

    if (coeff_file != "No File"){               /* Storing coefficient of each state */
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
            k++;

        }
        cfile.close();
    }
    else{
        for (int k = 0;k<states;k++){
            coeffs[k] = 1/pow(states,0.5);
        }
    }
    vector<vector<double>> wavefunction(n, vector<double>(n, 0.0));         /* Writing overall wave function to grid */
    ofstream file1;
    file1.open("outputs.txt", ios::app);
    file1 << endl;

    for (int i = 1; i < n - 1; i++) {
        file1 << endl;
        for (int k = 0; k < states; k++) {
            for (int j = 1; j < n - 1; j++) {

                file1 << abs(grid[i][j][k][0]) << ",";
                wavefunction[i][j] += abs(grid[i][j][k][0]);

            }
            file1 << "\t";

        }


        for (int j = 1; j < n - 1; j++) {
            file1 << wavefunction[i][j]/states << ",";
        }
    }
    file1.close();

}

void Ent_Grid::GetEntNextState(int n, double B, int states, int *signal) { /* Advancing grid by one time step */

    vector<vector<vector<complex<double>>>> neighbours(n, vector<vector<complex<double>>>(n, vector<complex<double>>(states, 0.0)));

    for (int k = 0; k < states; k++) {
        for (int i = 1; i < n - 1; i++) {                   // Checking number of neighbours
            for (int j = 1; j < n - 1; j++) {


                for (int a = -1; a < 2; a++) {              //Counting number of neighbours for each cell
                    for (int b = -1; b < 2; b++) {
                        if ((a == 0) && (b == 0)) {
                        } else {
                            neighbours[i][j][k] += grid[i + a][j + b][k][0];
                        }
                    }
                }
            }
        }
        for (int i = 1; i < n-1; i++) {                   // Checking number of neighbours and assigning new state
            for (int j = 1; j < n-1; j++) {
                int sign = 0;
                if (abs(grid[i][j][k][0]) != 0) {
                    grid[i][j][k][3].real(pow(real(grid[i][j][k][0]), 2) / norm(grid[i][j][k][0]));
                    grid[i][j][k][3].imag(pow(imag(grid[i][j][k][0]), 2) / norm(grid[i][j][k][0]));

                    if (grid[i][j][k][0].imag() < 0) {
                        sign = -1;
                    }
                    else{
                    }
                }
                if (abs(neighbours[i][j][k]) <= 1) {
                    grid[i][j][k][2].real(0); /* placeholder for alive state */
                    grid[i][j][k][1].real(1); /* dead state calculation */
                } else if ((abs(neighbours[i][j][k]) > 1) && (abs(neighbours[i][j][k]) <= 2)) {
                    grid[i][j][k][2].real((abs(neighbours[i][j][k])-1)*abs(grid[i][j][k][0]));
                    grid[i][j][k][1].real(B*(2-abs(neighbours[i][j][k]))*(real(grid[i][j][k][1])+abs(grid[i][j][k][0])) + (abs(neighbours[i][j][k])-1)*real(grid[i][j][k][1]));
                } else if ((abs(neighbours[i][j][k]) > 2) && (abs(neighbours[i][j][k]) <= 3)) {
                    grid[i][j][k][2].real(B*((3 - abs(neighbours[i][j][k])) * abs(grid[i][j][k][0])) + (abs(neighbours[i][j][k])-2)*(abs(grid[i][j][k][0]) + real(grid[i][j][k][1])));
                    grid[i][j][k][1].real(B*((3 - abs(neighbours[i][j][k])) * real(grid[i][j][k][1])));
                } else if ((abs(neighbours[i][j][k]) > 3) && (abs(neighbours[i][j][k]) <= 4)) {
                    grid[i][j][k][2].real(B*((4-abs(neighbours[i][j][k]))*(abs(grid[i][j][k][0]) + real(grid[i][j][k][1]))));
                    grid[i][j][k][1].real((abs(neighbours[i][j][k])-3)*(abs(grid[i][j][k][0]) + real(grid[i][j][k][1])));
                } else if ((abs(neighbours[i][j][k]) > 4)){
                    grid[i][j][k][2].real(0);
                    grid[i][j][k][1].real(1);
                }

                double norm_constant = real(grid[i][j][k][2]) + real(grid[i][j][k][1]);             /* Ensuring correct ratio of alive/dead coefficients */
                grid[i][j][k][2].real(grid[i][j][k][2].real()/norm_constant);
                grid[i][j][k][1].real(real(grid[i][j][k][1])/norm_constant);
                if (abs(grid[i][j][k][0]) != 0) {
                    grid[i][j][k][0].real(pow(real(grid[i][j][k][3]) * real(grid[i][j][k][2]), 0.5));
                    if(sign == -1) {
                        grid[i][j][k][0].imag(-pow(imag(grid[i][j][k][3]) * real(grid[i][j][k][2]), 0.5));
                    }
                    else{
                        grid[i][j][k][0].imag(pow(imag(grid[i][j][k][3]) * real(grid[i][j][k][2]), 0.5));
                    }
                }
                else{
                    if (neighbours[i][j][k].real() == 0){
                        if (neighbours[i][j][k].imag()<0) {
                            grid[i][j][k][0].imag(-grid[i][j][k][2].real());
                        }
                        else{
                            grid[i][j][k][0].imag(grid[i][j][k][2].real());
                        }
                    }
                    else if(neighbours[i][j][k].imag() == 0){
                        grid[i][j][k][0].real(grid[i][j][k][2].real());
                    }
                    else{
                        grid[i][j][k][0].real(grid[i][j][k][2].real()*pow(2,0.5)/2);

                        if (neighbours[i][j][k].imag()<0) {
                            grid[i][j][k][0].imag(-grid[i][j][k][0].real());
                        }
                        else{
                            grid[i][j][k][0].imag(grid[i][j][k][0].real());
                        }
                    }
                }

                grid[i][j][k][2].real(0);
            }


            }
        }
    vector<vector<double>> wavefunction(n, vector<double>(n, 0.0));
    ofstream myfile;
    myfile.open("outputs.txt", ios::app);
    myfile << endl;

    for (int i = 1; i < n - 1; i++) {                       /* Writing wave function to file */
        myfile << endl;
        for (int k = 0; k < states; k++) {
            for (int j = 1; j < n - 1; j++) {

                myfile << abs(grid[i][j][k][0]) << ",";
                wavefunction[i][j] += abs(grid[i][j][k][0]);

            }
            myfile << "\t";

        }


        for (int j = 1; j < n - 1; j++) {
            if (*signal == 0) {
                myfile << wavefunction[i][j] / states << ",";
            }
            else{
                myfile << wavefunction[i][j] << ",";
            }


        }
    }
    myfile.close();


    }



    double Ent_Grid::EntangleEntropy(int n, int states, string ent_file, int sub_cells) {  /* Function calculating entanglement entropy */

        vector<vector<double>> cellstates(states, vector<double>(pow(n,2),0));                                             /* State of each cell */
        vector<vector<int>> subspace_coords(sub_cells, vector<int>(2,0));                                                                /* Co-ords of subsspace */
        vector<int> cell_id(sub_cells,0);                                                                                                      /* Numbering each cell, and then denoting which ones are alive */
        vector<vector<int>> matched_pairs(states,vector<int>(states,1));                                                                 /* Used in trace calculation */
        vector<vector<complex<double>>> density_matrix(pow(2,pow(n,2)-sub_cells),                                     /* Density Matrix */
                                                       vector<complex<double>>(pow(2,pow(n,2)-sub_cells)));
        vector<complex<double>> density_vector(pow(pow(2,pow(n,2)-sub_cells),2),0);                  /* Density matrix as vector */
        vector<int> state_id(states,0);                                                                                                              /* Representation of state in binary string */
        complex<double> entropy;
        for (int k = 0; k < states; k++) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    cellstates[k][i*n+j] = grid[i+1][j+1][k][0].real();
                }
            }
        }
        ifstream myfile(ent_file);

        if (!myfile.is_open()) {
            cout << "No such entanglement file" << endl;
            return 0;
        }
        int z = 0;
        while (!myfile.eof()) {
            int i, j;
            myfile >> i >> j;
            subspace_coords[z][0] = i;
            subspace_coords[z][1] = j;
            z++;

        }
        for(int a = 0;a<sub_cells;a++){
            cell_id[a] = subspace_coords[a][0]*n +subspace_coords[a][1];
        }

        for (int k = 0;k<states;k++){
            for (int j = 0;j<states;j++){
                for (int a = 0;a<sub_cells;a++){
                    if(cellstates[k][cell_id[a]] == cellstates[j][cell_id[a]]){
                    }
                    else{
                        matched_pairs[k][j] = 0;
                    }

                }
            }
        }
        sort(cell_id.begin(),cell_id.end());
        for(int k = 0;k<states;k++){
            for(int a =0;a<sub_cells;a++) {
                cellstates[k].erase(cellstates[k].begin()+(cell_id[a]-a-1));        /* Erasing traced over cells */
            }
        }
        for (int k = 0;k<states;k++){
            for (int a=0;a<cellstates[0].size();a++){
            state_id[k] += cellstates[k][a]*pow(2,cellstates[0].size()-1-a);        /* Writing each state as binary string */
            }
        }

        int counter = 0;
        for(int j = 0; j < states; j++) {
            for(int k = 0;k<states;k++) {
                if(matched_pairs[j][k] == 1) {
                    density_matrix[pow(2,pow(n,2)-sub_cells)-1-state_id[k]][pow(2,pow(n,2)-sub_cells)-1-state_id[j]] += coeffs[k]*coeffs[j];
                }
            }
        }
        density_vector = density_matrix[0];
        for (int a = 1;a<=pow(2,pow(n,2)-sub_cells);a++){
            density_vector.insert(density_vector.end(),begin(density_matrix[a]),end(density_matrix[a]));
        }

        typedef Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> matrix;                                    /* Converting to matrix in Eigen library for entropy calculation */
        matrix dens_matrix = Eigen::Map<matrix>(density_vector.data(),pow(2,pow(n,2)-sub_cells),pow(2,pow(n,2)-sub_cells));

        cout << endl;
        VectorXcd e_values = dens_matrix.eigenvalues();
        for(int a =0;a<pow(2,pow(n,2)-sub_cells);a++){
            if (e_values(a) != 0.0){
                entropy += e_values(a)*(log(e_values(a)));
            }
        }

        cout << "Entropy of entanglement is" << -1.0*entropy << endl;
    return 0;
}












void Ent_Grid::Ent_Collapse(int n, int *signal, string coeff_file){
    double norm_factor = 0;
    int measured_state;
    for(int i = 0; i<coeffs.size();i++){
        norm_factor += abs(coeffs[i]);
    }

    if (coeff_file != "No File"){
        vector<int> ranges(coeffs.size());
        double ctotal;
        for (int m = 0; m<coeffs.size();m++){
            ctotal += (abs(coeffs[m])/norm_factor);
            ranges[m] = 10000*ctotal;
        }
        int r_num = (rand() % 10000);
        int index = 0;
        while(r_num>ranges[index]){
            index++;
            measured_state++;
        }

    }
    else {
        measured_state = (rand() % coeffs.size());
    }

    for (int l = 1; l < n - 1; l++) {
        for (int m = 0; m < coeffs.size(); m++) {
            for (int p = 1; p < n - 1; p++) {
                if (m != measured_state) {
                    grid[l][p][m][0].real(0);
                    grid[l][p][m][0].imag(0);
                    grid[l][p][m][1].real(1);
                    grid[l][p][m][1].imag(0);

                }
            }
        }
    }
    *signal = 1;
}


