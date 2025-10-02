#include <vector>
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

    //5-10 Iterate
    double p00prior=100.0;
    while (std::abs(P[0][0] - p00prior) > 0.1){
        p00prior=P[0][0];
        P = iterationP(P, intMatsR);


    }


    return 0;
}