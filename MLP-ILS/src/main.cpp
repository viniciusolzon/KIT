#include "MLP.h"
#include "readData.h"

int main(int argc, char *argv[]){

    double **matrizAdj;
    int dimension;

    readData(argc, argv, &dimension, &matrizAdj);
    
    Data d = {matrizAdj, dimension};
    MLP solution(d);
    
    solution.solve();
    solution.show_solution();

    return 0;
}
