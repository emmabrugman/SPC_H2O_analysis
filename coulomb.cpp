#include "coulomb.h"
#include <stdexcept>
#include <cmath>

double Coulomb::getPartialCharge(int Z) {
    if (Z == 6) return -0.82;  // Oxygen
    if (Z == 1) return +0.41;  // Hydrogen
    throw std::runtime_error("Unknown atomic number");
}

double Coulomb::computeEnergy(const arma::mat& positions, const std::vector<int>& A) {
    double energy = 0.0;

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = i + 1; j < A.size(); ++j) {
            double q1 = getPartialCharge(A[i]);
            double q2 = getPartialCharge(A[j]);
            double r = arma::norm(positions.col(i) - positions.col(j), 2);
            energy += COULOMB_CONSTANT * q1 * q2 / r;
        }
    }

    return energy;
}
