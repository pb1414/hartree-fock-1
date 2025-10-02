#include <vector>
#include <iostream>
#include <cmath> // Required for the pow() function
#include "params.h"
#include "integrals.h"


// =======================  Main  ==========================
void callHF(){
    //2-3 Calc matrices
    double R = 1.4;
    IntMats intMatsR = collect_integrals(R);
    //4 Guess P density matrix
    std::vector<std::vector<double>> P = {
        {0.0, 0.0},
        {0.0, 0.0}
    };
    

    //5-10 Call iterations

}


// std::vector<double> intgrl()


// =================== Diag fxn, MatMultfxn,  ====================

std::tuple<std::vector<std::vector<double>>,
                  std::vector<std::vector<double>>>
diag(const std::vector<std::vector<double>>& fmat) {
    double theta = PI / 4.0;

    if (std::abs(fmat[0][0] - fmat[1][1]) > 1e-10) {
        theta = 0.5 * std::atan(2.0 * fmat[0][1] / (fmat[0][0] - fmat[1][1]));
    }

    std::vector<std::vector<double>> cmat = {
        {std::cos(theta), std::sin(theta)},
        {std::sin(theta), -std::cos(theta)}
    };

    std::vector<std::vector<double>> emat(2, std::vector<double>(2, 0.0));

    emat[0][0] = fmat[0][0] * std::cos(theta) * std::cos(theta)
               + fmat[1][1] * std::sin(theta) * std::sin(theta)
               + fmat[0][1] * std::sin(2.0 * theta);

    emat[1][1] = fmat[1][1] * std::cos(theta) * std::cos(theta)
               + fmat[0][0] * std::sin(theta) * std::sin(theta)
               - fmat[0][1] * std::sin(2.0 * theta);

    // order eigvals correctly
    if (emat[1][1] > emat[0][0]) {
        std::swap(emat[0][0], emat[1][1]);
        std::swap(cmat[0][0], cmat[0][1]);
        std::swap(cmat[1][0], cmat[1][1]);
    }

    return {cmat, emat};
}
