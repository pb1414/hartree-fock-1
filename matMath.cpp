#include <vector>
#include <tuple>
#include <iostream>
#include <cmath> // Required for the pow() function
constexpr double PI = 3.14159265358979323;


// =================== Diag fxn ====================

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
    if (emat[1][1] < emat[0][0]) {
        std::swap(emat[0][0], emat[1][1]);
        std::swap(cmat[0][0], cmat[0][1]);
        std::swap(cmat[1][0], cmat[1][1]);
    }

    return {cmat, emat};
}


// =================== 2x2matMath ====================
std::vector<std::vector<double>> matMult(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B,
    int maxIndex) 
{
    // Initialize result matrix with zeros
    std::vector<std::vector<double>> C(2, std::vector<double>(2, 0.0));

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < maxIndex; k++) {
                C[i][j] += A[i][k] * B[k][j];
                if (maxIndex==1){
                    C[i][j] = C[i][j] * 2.0;  //This is activated when density matrix is being made
                }
            }
        }
    }
    return C;
}


