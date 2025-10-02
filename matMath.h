#pragma once
#include <vector>
#include <iostream>

std::tuple<std::vector<std::vector<double>>,
                  std::vector<std::vector<double>>>
diag(const std::vector<std::vector<double>>& fmat);

std::vector<std::vector<double>> matMult(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B,
    int maxIndex);