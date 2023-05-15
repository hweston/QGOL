#include <iostream>
#include <string>
#include "Grid.h"
#include "Ent_Grid.h"

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;



int main() {
  /*  clock_t start, end;                              For measuring time taken for code to run
    start = clock(); */



    int n = 60;                                     /* Gridsize */
    double B = (pow(2, 0.5) + 1);     /* Parameter controlling evolution rules */
    int steps;                                      /* Timesteps */
    int states = 4;                                 /* Size of grid vector 3rd dimension - corresponds to entangled states or alive/dead state in standard version*/
    int cells;                                      /* Stores total number of entangled cells */
    vector<int> collapsetimes(10,10000);      /* Time of measurements */
    int no_of_collapses = 0;                        /* No. of measurements */
    int signal = 0;                                 /* Measurement signal*/
    int sub_cells;                                  /* Size of subspace */

    string filename;                                /* Filename with initial seed */
    string coeff_file = "No File";                  /* Filename with coefficients for each state */
    string simulation;                              /* Specifies "entangled" or "standard simulation" */
    string collapse;                                /* Specifies whether any measurements take place */
    string ent;                                     /* Specifies whether entropy entanglement is calculated */
    string ent_file;                                /* File denoting coordiantes of cells in subspace */

    cout << "Type of Simulation: ";
    cin >> simulation;

    cout << "Do you want the wavefunction to collapse? (y/n) ";
    cin >> collapse;


    if (simulation == "standard") {

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
        if (collapse == "y") {
            cout << "After which timestep?: ";
            cin >> collapsetimes[0];
        }
        string coeff;
        cout << "Are coefficients for each state being chosen? (y/n): ";
        cin >> coeff;
        if(coeff == "y"){
            coeff_file = "temp";
        }

        if (coeff_file != "No File") {
            cout << "Please enter filename listing coeffient values: ";
            cin >> coeff_file;
        }

        cout << "How many states are entangled?: ";
        cin >> states;
        cout << "How many total cells are entangled?: ";
        cin >> cells;
    }

    vector<int> cellnum(states);                 /* Stores number of cells in each entangled state */

    if(simulation == "entangled"){                  /* Checking and storing number of cells entangled*/
        int checkcells = 0;
        for (int i = 0; i < states; i++) {
                cout << "How many cells are in state " << i + 1
                     << "?: ";
                cin >> cellnum[i];
                checkcells += cellnum[i];
                if (i == states - 1) {
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

    if (simulation == "standard") {                             /* Performs standard (non-entangled) simulation */

        Grid state = Grid(filename, n + 2);

        for (int i = 0; i < steps; i++) {

            for (int a = 0; a < no_of_collapses; a++)
                if (i == collapsetimes[a]) {
                    for (int j = 0; j < n; j++) {
                        for (int k = 0; k < n; k++) {
                            state.Collapse(j+1, k+1);
                        }
                    }
                }
            state.GetNextState(n + 2, B);
        }
    }

    if (simulation == "entangled") {                                                        /* Performs entangled simulation */
        Ent_Grid state = Ent_Grid(filename, coeff_file, cells, cellnum, n + 2, states);
        for (int i = 0; i < steps; i++) {
            if (i == collapsetimes[0]) {
                state.Ent_Collapse(n + 2, &signal, coeff_file);
            }
            state.GetEntNextState(n + 2, B, states, &signal);
        }
        cout << "Calculate Entanglement Entropy? (y/n): ";
        cin >> ent;
        if (ent == "y"){
            cout << "Please enter filename which cell coordinates of subspace traced over: ";
            cin >> ent_file;
            cout << "How many cells in the subspace being traced over?: ";
            cin >> sub_cells;
            state.EntangleEntropy(n, states, ent_file, sub_cells);
        }

    }
    /*end = clock();                                                Code for measuring speed of program
    double tt = double(end - start)/double(CLOCKS_PER_SEC);
    cout << fixed << "Time taken to run: " << tt;
    return 0;
     */
}
