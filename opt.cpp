#include <vector>
#include <iostream>
#include <cmath> // Required for the pow() function
constexpr double PI=3.14159265358979323;

// Fxn headers. Defined below
double Sfxn(double alpha,double beta, double r); double T(double alpha,double beta, double r); double V(double alpha, double beta, double rab2, double rcp2, double zc);
double twoe(double alpha, double beta, double cgamma, double ddelta, double rab2, double rcd2, double rpq2);


// Setup for HF matrix elements
constexpr int nBas = 2;   //Code is very simple and assumes this, also same as number of nuclei
constexpr int nPrim = 3;
constexpr double lam = 1.6875;
constexpr double R = 1.4632;
const double alphas[nBas][nPrim] = {
    {0.168856, 0.623913, 3.42525},
    {0.168856 * lam * lam, 0.623913 * lam * lam, 3.42525 * lam * lam}
};
const double coefs[nPrim] = {0.444635,0.535328,0.154329};
const double Zs[nBas] = {1.0 , 1.0};

double main(){
    //2 Calc matrices
    
    //3 Get X

    //4 Guess P

    //5-10 Call iterations

}

// std::vector<double> intgrl()

// ===================  Collecting integrals ===================

void collect_integrals(){
    //Generate S
    std::vector<std::vector<double>> S= {
        {1.0,0.0},
        {0.0,1.0}
    };
    for (int i = 0; i<nPrim; i++){
        for (int j = 0; j<nPrim; j++){
            S[1][0]=S[1][0] + Sfxn(alphas[0][i],alphas[1][j], R);
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
                    double rab2 = (bi == bj) ? 0.0 : R*R;

                    Hcore[bi][bj] += T(alpha, beta, rab2);

                    for (int c = 0; c < nBas; c++) {
                        double Rc = (c == 0) ? 0.0 : R;
                        double P = (alpha*((bi==0)?0.0:R) + beta*((bj==0)?0.0:R)) / (alpha+beta);
                        double rcp2 = (P - Rc)*(P - Rc);
                        Hcore[bi][bj] += V(alpha, beta, rab2, rcp2, Zs[c]);
                    }

                    for (int bk = 0; bk < nBas; bk++) {
                        for (int bl = 0; bl < nBas; bl++) {
                            for (int k = 0; k < nPrim; k++) {
                                for (int l = 0; l < nPrim; l++) {
                                    double cgamma = alphas[bk][k];
                                    double ddelta = alphas[bl][l];
                                    double rcd2 = (bk == bl) ? 0.0 : R*R;
                                    double P = (alpha*((bi==0)?0.0:R) + beta*((bj==0)?0.0:R)) / (alpha+beta);
                                    double Q = (cgamma*((bk==0)?0.0:R) + ddelta*((bl==0)?0.0:R)) / (cgamma+ddelta);
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


    //
}


// ================= Basis fxn integrals ===================

// Overlap integral
double Sfxn(double alpha,double beta, double r){
    double ovlap= PI/(pow((alpha+beta),1.5)) * std::exp(-alpha*beta*r*r/(alpha+beta));
    return ovlap;
}

// KE integral
double T(double alpha,double beta, double r){
    double res= PI/(pow((alpha+beta),1.5)) * std::exp(-alpha*beta*R*R/(alpha+beta));
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


// =================== Diagonalizing fxn ====================

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
