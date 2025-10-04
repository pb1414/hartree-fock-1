#include <vector>
#include <iostream>
#include <tuple>
#include <cmath>
#include "params.h"
#include "integrals.h"
#include "iteration.h"
#include "matMath.h"

// ==================== Guessing P ====================

std::vector<std::vector<double>> guessP(){
    std::vector<std::vector<double>> guess_naive = {
        {0.0, 0.0},
        {0.0, 0.0}
    };
    return guess_naive;
}

// ==================== Single scf iteration ================

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> iterationP(std::vector<std::vector<double>> P_in, IntMats intMatsR){
    //5-6 Add G terms to Hcore to get F
    std::vector<std::vector<double>> F= intMatsR.Hcore;
    for (int bi = 0; bi < nBas; bi++) {
        for (int bj = 0; bj < nBas; bj++) {
            for (int bk = 0; bk < nBas; bk++) {
                for (int bl = 0; bl < nBas; bl++) {
                    F[bi][bj] = F[bi][bj] + P_in[bk][bl]*(intMatsR.TwoEIntegrals[bi][bj][bl][bk] - 0.5 * intMatsR.TwoEIntegrals[bi][bk][bl][bj]);
                }
            }
        }
    }


    //7 Get F'
    std::vector<std::vector<double>> mat1= matMult(intMatsR.XT, F, 2);
    std::vector<std::vector<double>> Fprime= matMult(mat1, intMatsR.X, 2);

    //8 Diagonalize F' to get the expansions of orbitals
    std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> Cs_Es  = diag(Fprime);

    //9 Get C back in correct basis
    std::vector<std::vector<double>> C = matMult(intMatsR.X, std::get<0>(Cs_Es), 2);

    //10 Get density matrix
        std::vector<std::vector<double>> CT = {
            {C[0][0], C[1][0]},
            {C[0][1], C[1][1]}
        };
    std::vector<std::vector<double>> P_new = matMult(C, CT, 1);

    return {P_new, std::get<1>(Cs_Es)};
}