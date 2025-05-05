#pragma once

#include "spc.h"
#include "forces.h"
#include <armadillo>

class CalculateHessian {
private:
    const SPC& spc_potential; // Reference to SPC potential
    int num_atoms;           // Number of atoms
    const CalculateForce& force_calculator; // Reference to force calculator (not used anymore)

public:
    // Constructor
    CalculateHessian(const SPC& spc, const CalculateForce& force_calc, int atoms);

    // Calculate numerical Hessian using central difference on forces
    void calculateCentralDifferenceHessian(double h, arma::mat& hessian) const;
};