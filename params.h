#pragma once
#include <vector>

// Setup for HF matrix elements
constexpr double PI = 3.14159265358979323;
constexpr int nBas = 2;
constexpr int nPrim = 3;
constexpr double lam = 1.6875;

const double alphas[nBas][nPrim] = {
    {0.168856 * lam * lam, 0.623913 * lam * lam, 3.42525 * lam * lam},
    {0.168856, 0.623913, 3.42525}
};

const double coefs[nPrim] = {0.444635, 0.535328, 0.154329};

const double coefsNorm[nBas][nPrim] = {
    {std::pow(2.0 * alphas[0][0] / PI, 0.75)*coefs[0],
    std::pow(2.0 * alphas[0][1] / PI, 0.75)*coefs[1],
    std::pow(2.0 * alphas[0][2] / PI, 0.75)*coefs[2]},
    {std::pow(2.0 * alphas[1][0] / PI, 0.75)*coefs[0],
    std::pow(2.0 * alphas[1][1] / PI, 0.75)*coefs[1],
    std::pow(2.0 * alphas[1][2] / PI, 0.75)*coefs[2]}
};
const double Zs[nBas] = {2.0, 1.0};
