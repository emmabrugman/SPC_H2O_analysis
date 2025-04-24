#include "forces.h"
#include <iostream>
#include <iomanip>

// Constructor
CalculateForce::CalculateForce(const SPC& spc, int atoms)
    : spcPotential(spc), num_atoms(atoms) {}

// Compute energy for a given displacement
double CalculateForce::calculateDisplacementEnergy(const arma::mat& newPositions) const {
    SPC spcDisplaced(newPositions, spcPotential.getAtomicNumbers());
    return spcDisplaced.getTotalEnergy();
}

// Central Difference: f'(x) ≈ -(f(x+h) - f(x-h)) / 2h
double CalculateForce::calculateCentralDifferences(double forward_energy, double backward_energy, double h) const {
    return -(forward_energy - backward_energy) / (2 * h);
}

// Compute central difference forces for each atom and coordinate
void CalculateForce::calculateCentralDifferenceForces(double h, arma::mat& central_forces) const {
    central_forces.zeros(3, num_atoms);

    const arma::mat& initialPositions = spcPotential.getPositions();
    double original_energy = spcPotential.getTotalEnergy();

    // Loop over each atom and each coordinate
    for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
            arma::mat forward_positions = initialPositions;
            arma::mat backward_positions = initialPositions;

            // Displace in the positive direction
            forward_positions(coordinate_index, atom_idx) += h;

            // Displace in the negative direction
            backward_positions(coordinate_index, atom_idx) -= h;

            double forward_energy = calculateDisplacementEnergy(forward_positions);
            double backward_energy = calculateDisplacementEnergy(backward_positions);

            // Compute force using central difference approximation
            central_forces(coordinate_index, atom_idx) = calculateCentralDifferences(forward_energy, backward_energy, h);
        }
    }
}