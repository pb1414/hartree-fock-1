#include <vector>
#include <iostream>
#include <cmath> // Required for the pow() function
constexpr double PI=3.14159265358979323;
std::vector<double> intgrl();


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
const double coefs[nBas] = {0.444635,0.535328,0.154329};

double main(){
    //2 Calc integrals 
    
    //3 Get X

    //4 Guess P

    //5-10 Call iterations

}


// std::vector<double> intgrl()

void collect_ints();




// ================= Basis fxn integrals ===================

// Overlap integral
double S(double alpha,double beta, double r){
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

inline std::tuple<std::vector<std::vector<double>>,
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
