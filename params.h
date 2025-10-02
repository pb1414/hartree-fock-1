#pragma once
#include <vector>

// Setup for HF matrix elements
constexpr double PI = 3.14159265358979323;
constexpr int nBas = 2;
constexpr int nPrim = 3;
constexpr double lam = 1.6875;

const double alphas[nBas][nPrim] = {
    {0.168856, 0.623913, 3.42525},
    {0.168856 * lam * lam, 0.623913 * lam * lam, 3.42525 * lam * lam}
};

const double coefs[nPrim] = {0.444635, 0.535328, 0.154329};
const double Zs[nBas] = {1.0, 1.0};
