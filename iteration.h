#pragma once
#include <vector>
#include <iostream>
#include "integrals.h"

std::vector<std::vector<double>> guessP();
std::vector<std::vector<double>> iterationP(std::vector<std::vector<double>> P_in, IntMats intMatsR);
