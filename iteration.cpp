#include <vector>
#include <iostream>
#include <cmath>
#include "params.h"
#include "integrals.h"
#include "iteration.h"
#include "matMath.h"

// ==================== Guessing P ====================

std::vector<std::vector<double>> guessP(){
    std::vector<std::vector<double>> guess_naive = {
        {0.25, 0,25},
        {0.25,0.25}
    };
    return guess_naive;
}

// ==================== Single scf iteration ================

std::vector<std::vector<double>> iterationP(std::vector<std::vector<double>> P_in, IntMats intMatsR){

    //5-6 Add G terms to Hcore to get F
    std::vector<std::vector<double>> F= intMatsR.Hcore;
    for (int bi = 0; bi < nBas; bi++) {
        for (int bj = 0; bj < nBas; bj++) {
            for (int bk = 0; bk < nBas; bi++) {
                for (int bl = 0; bl < nBas; bj++) {
                    F[bi][bj] = F[bi][bj] + P_in[bi][bj]*intMatsR.TwoEIntegrals[bi][bj][bk][bl];
                }
            }
        }
    }

    //7 Get F'
    std::vector<std::vector<double>> mat1= matMult(intMatsR.XT,F);
    std::vector<std::vector<double>> Fprime= matMult(mat1,intMatsR.X);

    //8 Diagonalize F' to get the expansions of orbitals
    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> Cs_Es  = diag(Fprime);

    //9 Get C back in correct basis
    std::vector<std::vector<double>> C = matMult(intMatsR.X, std::get<0>(Cs_Es));



    std::vector<std::vector<double>> P_new = {
        {0.25, 0,25},
        {0.25,0.25}
    };

    return P_new;
}