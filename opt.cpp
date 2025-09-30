#include <vector>
#include <iostream>


// Setup for HF matrix elements
constexpr int nEl = 2;
constexpr int nBas = 3;
constexpr double lam = 1.6875;
constexpr double R = 1.4632;
const double alphas1[nBas] = {0.168856, 0.623913, 3.42525};
const double alphas2[nBas] = {
    alphas1[0] * lam * lam,
    alphas1[1] * lam * lam,
    alphas1[2] * lam * lam
};


// for (int i = 0; i < nBas; i++) {
//     alphas2[i] = alphas1[i] * lam * lam;
// }




// std::vector<double> PreCoeffs(std::vector<double>(nBasis,0.0))

// int
// std::vector<std::vector<double>> alphas(nBas, std::vector<double>(nBas, 0.0));
// alphas[0][0] = 0.13;



double S(double r1, double r2){
    double combo = r1+r2;
    return combo;
}