#include <vector>
#include <tuple>
#include <iostream>
#include <cmath> // Required for the pow() function
#include "params.h"
#include "integrals.h"
#include "iteration.h"


// =======================  Main  ==========================
int main(){
    //2-3 Calc matrices
    double R = 1.4632;
    IntMats intMatsR = collect_integrals(R);

    //4 Guess P density matrix
    std::vector<std::vector<double>> P = guessP();
    std::vector<std::vector<double>> P2 = P;

    std::cout << "Started" << std::endl;

    //5-10 Iterate
    std::vector<std::vector<double>> es(2, std::vector<double>(2, 0.0));
    double p00prior=1.000;
    
    while (abs(P[0][0]-p00prior)>0.001){
        p00prior=P[0][0];
        std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> PandEs = iterationP(P, intMatsR);
        P= std::get<0>(PandEs);
        es= std::get<1>(PandEs);
        std::cout << "New eigenvalues: " << es[0][0] << ", " << es[1][1] << std::endl;
    }


    //12 Report results
    
    std::cout << "Finished. Final eigenvalue matrix:" << std::endl;
    for (const auto& row : es) {
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}