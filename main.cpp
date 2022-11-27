#include <iostream>
#include "Grid.h"
#include "GetState.h"

using namespace std;


int main()
{
    string filename;
    int n = 25;
    int steps;
    cout << "Enter filename with initial grid: ";
    cin >> filename;
    cout << "Enter number of steps of simulation: ";
    cin >> steps;
    Grid state = Grid(filename, n + 2);

    for(int i =0;i<steps;i++){
        state.GetNextState(n+2);
    }
}





