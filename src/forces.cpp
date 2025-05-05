#include "forces.h"
#include <iostream>
#include <iomanip>

// Constructor
CalculateForce::CalculateForce(const SPC& spc, int atoms)
    : spc_potential(spc), num_atoms(atoms) {}

// Compute energy for a given displacement
double CalculateForce::calculateDisplacementEnergy(const arma::mat& new_positions) const {
    SPC spc_displaced(new_positions, spc_potential.getAtomicNumbers());
    return spc_displaced.getTotalEnergy();
}

// Central Difference: f'(x) ≈ -(f(x+h) - f(x-h)) / 2h
double CalculateForce::calculateCentralDifferences(double forward_energy, double backward_energy, double h) const {
    return -(forward_energy - backward_energy) / (2 * h);
}

// Compute central difference forces for each atom and coordinate
void CalculateForce::calculateCentralDifferenceForces(double h, arma::mat& central_forces) const {
    central_forces.zeros(3, num_atoms);

    const arma::mat& initial_positions = spc_potential.getPositions();
    double original_energy = spc_potential.getTotalEnergy();

    // Loop over each atom and each coordinate
    for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        for (int coord_idx = 0; coord_idx < 3; ++coord_idx) {
            arma::mat forward_positions = initial_positions;
            arma::mat backward_positions = initial_positions;

            // Displace in the positive direction
            forward_positions(coord_idx, atom_idx) += h;

            // Displace in the negative direction
            backward_positions(coord_idx, atom_idx) -= h;

            double forward_energy = calculateDisplacementEnergy(forward_positions);
            double backward_energy = calculateDisplacementEnergy(backward_positions);

            // Compute force using central difference approximation
            central_forces(coord_idx, atom_idx) = calculateCentralDifferences(forward_energy, backward_energy, h);
        }
    }
}