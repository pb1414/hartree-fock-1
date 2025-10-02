#include <vector>
#include <iostream>
#include <cmath> // Required for the pow() function
#include "params.h"
double Sfxn(double alpha,double beta, double r); double T(double alpha,double beta, double r); double V(double alpha, double beta, double rab2, double rcp2, double zc);
double twoe(double alpha, double beta, double cgamma, double ddelta, double rab2, double rcd2, double rpq2); struct IntMats;


// ===================  Collecting integrals ===================

IntMats collect_integrals(double r){
    //Generate S
    std::vector<std::vector<double>> S= {
        {1.0,0.0},
        {0.0,1.0}
    };
    for (int i = 0; i<nPrim; i++){
        for (int j = 0; j<nPrim; j++){
            S[1][0]=S[1][0] + Sfxn(alphas[0][i],alphas[1][j], r);
        }
    }
    S[0][1]=S[1][0];

    // Generate X
    double s12 = S[0][1];  // assuming S is 2x2

    std::vector<std::vector<double>> X = {
        {1.0 / std::sqrt(2.0 * (1.0 + s12)),  1.0 / std::sqrt(2.0 * (1.0 - s12))},
        {1.0 / std::sqrt(2.0 * (1.0 + s12)), -1.0 / std::sqrt(2.0 * (1.0 - s12))}
    };

    std::vector<std::vector<double>> XT = {
        {X[0][0], X[1][0]},
        {X[0][1], X[1][1]}
    };

    // Generate Hcore and the 2e matrix,  bi and bj go over bases
    std::vector<std::vector<double>> Hcore(nBas, std::vector<double>(nBas, 0.0));
    std::vector<std::vector<std::vector<std::vector<double>>>> TwoEIntegrals(
        nBas, std::vector<std::vector<std::vector<double>>>(
                nBas, std::vector<std::vector<double>>(
                            nBas, std::vector<double>(nBas, 0.0))));

    for (int bi = 0; bi < nBas; bi++) {
        for (int bj = 0; bj < nBas; bj++) {
            for (int i = 0; i < nPrim; i++) {
                for (int j = 0; j < nPrim; j++) {
                    double alpha = alphas[bi][i];
                    double beta  = alphas[bj][j];
                    double rab2 = (bi == bj) ? 0.0 : r*r;

                    Hcore[bi][bj] += T(alpha, beta, rab2);

                    for (int c = 0; c < nBas; c++) {
                        double Rc = (c == 0) ? 0.0 : r;
                        double P = (alpha*((bi==0)?0.0:r) + beta*((bj==0)?0.0:r)) / (alpha+beta);
                        double rcp2 = (P - Rc)*(P - Rc);
                        Hcore[bi][bj] += V(alpha, beta, rab2, rcp2, Zs[c]);
                    }

                    for (int bk = 0; bk < nBas; bk++) {
                        for (int bl = 0; bl < nBas; bl++) {
                            for (int k = 0; k < nPrim; k++) {
                                for (int l = 0; l < nPrim; l++) {
                                    double cgamma = alphas[bk][k];
                                    double ddelta = alphas[bl][l];
                                    double rcd2 = (bk == bl) ? 0.0 : r*r;
                                    double P = (alpha*((bi==0)?0.0:r) + beta*((bj==0)?0.0:r)) / (alpha+beta);
                                    double Q = (cgamma*((bk==0)?0.0:r) + ddelta*((bl==0)?0.0:r)) / (cgamma+ddelta);
                                    double rpq2 = (P - Q)*(P - Q);
                                    TwoEIntegrals[bi][bj][bk][bl] += twoe(alpha, beta, cgamma, ddelta, rab2, rcd2, rpq2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return {S, X, XT, Hcore, TwoEIntegrals};
}


// ================= Basis fxn integrals ===================

// Overlap integral
double Sfxn(double alpha,double beta, double r){
    double ovlap= PI/(pow((alpha+beta),1.5)) * std::exp(-alpha*beta*r*r/(alpha+beta));
    return ovlap;
}

// KE integral
double T(double alpha,double beta, double r){
    double res= PI/(pow((alpha+beta),1.5)) * std::exp(-alpha*beta*r*r/(alpha+beta));
    return res;
}

// Boys helper
double boys(double t) {
    if (t < 1e-8) {
        // limit F0(t → 0) = 1
        return 1.0;
    } else {
        return 0.5 * std::sqrt(PI / t) * std::erf(std::sqrt(t));
    }
}

// Nuc attraction integral
double V(double alpha, double beta, double rab2, double rcp2, double zc){
    double res = -2.0/zc*PI/(alpha+beta)*boys((alpha+beta)*rcp2)*std::exp(-1.0*alpha*beta*rab2/(alpha+beta));
    return res;
}

// 2e integral
double twoe(double alpha, double beta, double cgamma, double ddelta,
                   double rab2, double rcd2, double rpq2) {
    double denom = (alpha + beta) * (cgamma + ddelta) * std::sqrt(alpha + beta + cgamma + ddelta);
    double pre = 2.0 * std::pow(PI, 2.5) / denom;
    double t = (alpha + beta) * (cgamma + ddelta) * rpq2 / (alpha + beta + cgamma + ddelta);
    double expo = std::exp(-alpha * beta * rab2 / (alpha + beta)
                           - cgamma * ddelta * rcd2 / (cgamma + ddelta));
    return pre * boys(t) * expo;
}


// ================= IntMats Struct ===================


struct IntMats {
    std::vector<std::vector<double>> S;
    std::vector<std::vector<double>> X;
    std::vector<std::vector<double>> XT;
    std::vector<std::vector<double>> Hcore;
    std::vector<std::vector<std::vector<std::vector<double>>>> TwoEIntegrals;
};