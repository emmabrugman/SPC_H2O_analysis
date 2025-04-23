#pragma once

#include <armadillo>
#include <vector>

class Coulomb {
public:
    static constexpr double COULOMB_CONSTANT = 332.0637; // kJ·Å/(mol·e²)

    // Compute total Coulomb energy
    static double computeEnergy(const arma::mat & positions, const std::vector<int> & atomic_numbers);

    // SPC partial charge based on atomic number
    static double getPartialCharge(int atomic_number);
};
