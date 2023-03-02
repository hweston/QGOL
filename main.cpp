#include <iostream>
#include <string>
#include "Grid.h"


using namespace std;


int main()
{
    string filename;
    int n = 6;
    double B = 10;
    int steps;
    int doe = 1;
    string simulation;
    string entangled = "entangled";
    string standard = "standard";

    cout << "Type of Simulation: ";
    cin >> simulation;
    if (simulation == entangled){
        cout << "Enter number of entangled cells: ";
        cin >> doe;
    }
    cout << "Enter filename with initial grid: ";
    cin >> filename;
    cout << "Enter number of steps of simulation: ";
    cin >> steps;

    Grid state = Grid(filename, n + 2, B,simulation,doe);

    if (simulation == standard) {

                for (int i =0;i<steps;i++){
                    state.GetNextState(n + 2, B, simulation,1);
                }
        }
    if (simulation == entangled) {
        for (int i =0;i<steps;i++){
            state.GetNextState(n + 2, B,simulation, doe);
        }
    }


    else {
            cout << "Error" << endl;
        }
        }







