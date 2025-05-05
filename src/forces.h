#pragma once

#include "spc.h"
#include <armadillo>

class CalculateForce {
private:
    const SPC& spc_potential; // Reference to SPC potential
    int num_atoms;           // Number of atoms

    // Compute energy for a given displacement
    double calculateDisplacementEnergy(const arma::mat& new_positions) const;

    // Central difference approximation
    double calculateCentralDifferences(double forward_energy, double backward_energy, double h) const;

public:
    // Constructor
    CalculateForce(const SPC& spc, int atoms);

    // Calculate finite difference forces (central difference only)
    void calculateCentralDifferenceForces(double h, arma::mat& central_forces) const;
};