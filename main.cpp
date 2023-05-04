#include <iostream>
#include <string>
#include "Grid.h"
#include <Eigen/Dense>



using namespace std;
using namespace Eigen;



int main() {

    int n = 30;                     /* Gridsize */
    double B = (pow(2,0.5)+1);                 /* Parameter controlling evolution rules */
    int steps;                      /* Timesteps */
    int z = 4;                      /* Size of grid vector 3rd dimension - corresponds to entangled states/alive/dead state */
    int signal = 0;                 /* Signals wavefunction collapse - NB different for entangled vs. standard */
    int cells;                      /* Stores total number of cells */
    vector<int> collapsetimes(10);
    int no_of_collapses;

    string filename;                /* Filename with initial seed */
    string coeff_file;              /* Filename with coefficients for each state */
    string simulation;              /* Specifies "entangled" or "standard simulation" */
    string etype = "NA";            /* All entanglement or specified number of states */
    string coeff;                   /* Specifies whether coefficients are chosen or not */
    string collapse;

    for (int i = 0; i < collapsetimes.size(); i++) {
        collapsetimes[i] = 10000;
    }

    cout << "Type of Simulation: ";
    cin >> simulation;
    if (simulation == "standard") {
        cout << "Do you want the wavefunction to collapse? (y/n) ";
        cin >> collapse;
        if (collapse == "y") {
            cout << "How many times?: ";
            cin >> no_of_collapses;
            for (int i = 0; i < no_of_collapses; i++) {
                cout << "Time of collapse " << i + 1 << ": ";
                cin >> collapsetimes[i];
            }

        }
    }
    if (simulation == "entangled") {
        cout << "All states(y) or a chosen number of states(n)?: ";
        cin >> etype;
        cout << "Are coefficients for each state being chosen? (y/n): ";
        cin >> coeff;
        if (coeff == "y") {
            cout << "Please enter filename listing coeffient values: ";
            cin >> coeff_file;
        }
        if (etype == "n") {
            cout << "How many states are entangled?: ";
            cin >> z;
        } else {

            cout << "Enter number of entangled states: ";
            cin >> z;
        }

    }
    vector<int> cellnum(z);   /* Stores number of cells in each entangled state */
    if (etype == "n") {
        cout << "How many total cells are entangled?: ";
        cin >> cells;
        int checkcells = 0;
        for (int i = 0; i < z; i++) {
            cout << "How many cells are in state " << i + 1
                 << "?: ";
            cin >> cellnum[i];
            checkcells += cellnum[i];
            if (i == z - 1) {
                if (checkcells != cells) {
                    cout
                            << "Error: Your total cells doesn't match the sum of cells from each state. Please enter cell numbers again: ";
                    i = 0;
                }
            }
        }
    }

    cout << "Enter filename with initial grid: ";
    cin >> filename;
    cout << "Enter number of steps of simulation: ";
    cin >> steps;


    Grid state = Grid(filename, coeff_file, coeff, cells, cellnum, n + 2, simulation, etype, z);

    if (simulation == "standard") {
        for (int i = 0; i < steps; i++) {
            for (int a = 0; a < no_of_collapses; a++)
                if (i == collapsetimes[a]) {
                    for (int j = 1; j < n + 1; j++) {
                        for (int k = 1; k < n + 1; k++) {
                            state.Collapse(j, k);
                        }
                    }
                }
            state.GetNextState(n + 2, B, simulation, coeff, 1, etype, &signal);
        }
    }
    if (simulation == "entangled") {

        for (int i = 0; i < steps; i++) {
            state.GetNextState(n + 2, B, simulation, coeff, z, etype, &signal);
        }
    }

/*
    state.EntangEntropy(z,cellnum,cells);
*/
}

/*
 * Change Entang entropy function.
 * Form of wavefunction =
 * psi = Aphi1 + Bphi2
 * - need cells vector to be 2^(number of entangled cells) and then write each entangled state in binary to give state (phi1)
 * - then same for phi 2.
 * - NOT one state or the other (although could be structured this way - don't think like this)
 * Will probably struggle with more than 12 (?) entangled cells, this can be tested.
 * try and work out how to calculate reduced density matrix.
 * Probably best to write "cells" vector so that partitioned areas A and B are not mixed.
 *
 */




