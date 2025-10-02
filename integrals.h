#pragma once
#include <vector>
#include <iostream>

//Integral grouping
struct IntMats {
    std::vector<std::vector<double>> S;
    std::vector<std::vector<double>> X;
    std::vector<std::vector<double>> XT;
    std::vector<std::vector<double>> Hcore;
    std::vector<std::vector<std::vector<std::vector<double>>>> TwoEIntegrals;
};
IntMats collect_integrals(double r);
double Sfxn(double alpha, double beta, double r);
double T(double alpha, double beta, double r);
double boys(double t);
double V(double alpha, double beta, double rab2, double rcp2, double zc);
double twoe(double alpha, double beta, double cgamma, double ddelta,
            double rab2, double rcd2, double rpq2);