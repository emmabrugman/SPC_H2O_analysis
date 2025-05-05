#include "hessian.h"
#include <iostream>
#include <iomanip>

// Constructor
CalculateHessian::CalculateHessian(const SPC& spc, const CalculateForce& force_calc, int atoms)
    : spc_potential(spc), force_calculator(force_calc), num_atoms(atoms) {}

// The Hessian H_ij = ∂²E/(∂x_i ∂x_j) = -∂F_i/∂x_j
// This is calculated as: H_ij ≈ -(F_i(x_j+h) - F_i(x_j-h))/(2h)
void CalculateHessian::calculateCentralDifferenceHessian(double h, arma::mat& hessian) const {
    // Hessian is (3*N) x (3*N) where N is number of atoms
    // Each atom has 3 coordinates (x, y, z)
    hessian.zeros(3 * num_atoms, 3 * num_atoms);
    
    // Pre-allocate matrices to store forces at displaced positions
    // Each is 3 x N (3 coordinates per atom)
    const arma::mat& positions = spc_potential.getPositions();
    arma::mat forward_forces(3, num_atoms);
    arma::mat backward_forces(3, num_atoms);
    

    // Loop through each degree of freedom 
    // Each atom has 3 (x, y, z coordinates)
    for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        for (int coord_idx = 0; coord_idx < 3; ++coord_idx) {

            int row_idx = 3 * atom_idx + coord_idx; // Row index in Hessian matrix
            
            // FORWARD DISPLACEMENT: displace current coordinate by +h
            arma::mat forward_positions = positions;  // Copy original positions
            forward_positions(coord_idx, atom_idx) += h;  // Displace by +h
            SPC spc_forward(forward_positions, spc_potential.getAtomicNumbers()); // / Create SPC object with displaced positions
            CalculateForce force_calc_forward(spc_forward, num_atoms); // Temporary force calculator for displaced system
            force_calc_forward.calculateCentralDifferenceForces(h, forward_forces); // / Calculate forces at +h displaced position
            
            // BACKWARD DISPLACEMENT: displace current coordinate by -h
            arma::mat backward_positions = positions;
            backward_positions(coord_idx, atom_idx) -= h;
            SPC spc_backward(backward_positions, spc_potential.getAtomicNumbers());
            CalculateForce force_calc_backward(spc_backward, num_atoms);
            force_calc_backward.calculateCentralDifferenceForces(h, backward_forces);
            
            // COMPUTE HESSIAN ELEMENTS
            // Convert 3 x N force matrices to flat vectors of size 3N
            arma::vec forward_flat = arma::vectorise(forward_forces);
            arma::vec backward_flat = arma::vectorise(backward_forces);
            
            // Apply central difference formula to get Hessian column
            // Negative sign accounts for H = -∂F/∂x
            arma::vec hessian_column = -(forward_flat - backward_flat) / (2 * h);
            
            hessian.col(row_idx) = hessian_column;
            hessian.row(row_idx) = hessian_column.t();  // Symmetry
        }
    }
}