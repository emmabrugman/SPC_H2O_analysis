#include "lennard_jones.h"
#include <stdexcept>
#include <cmath>

// LJ constants
const double LennardJones::SIGMA_AU = 3.166; // Å
const double LennardJones::EPSILON_AU = 0.650; // kJ/mol

// Constructor: initialize positions and atomic numbers
LennardJones::LennardJones(const arma::mat& pos, const std::vector<int>& a)
    : positions(pos), A(a) {
}

// Compute Euclidean distance between two atoms
double LennardJones::calculateDistance(int first_atom_index, int second_atom_index) const {
    arma::vec diff = positions.col(first_atom_index) - positions.col(second_atom_index);
    return arma::norm(diff, 2);
}

double LennardJones::calculateR6Term(double r) const {
    return std::pow(SIGMA_AU / r, 6);
}

double LennardJones::calculateR12Term(double r6_term) const {
    return r6_term * r6_term;
}

double LennardJones::calculateLJ(double r_ij) const {
    double r6_term = calculateR6Term(r_ij);
    double r12_term = calculateR12Term(r6_term);
    return EPSILON_AU * (r12_term - 2 * r6_term);
}

double LennardJones::getLJEnergy() const {
    double lj_energy = 0.0;
    for (int atom1_idx = 0; atom1_idx < A.size(); ++atom1_idx) {
        for (int atom2_idx = atom1_idx + 1; atom2_idx < A.size(); ++atom2_idx) {
            // Apply LJ only if both atoms are oxygen (Z=6)
            if (A[atom1_idx] == 6 && A[atom2_idx] == 6) {
                double r_ij = calculateDistance(atom1_idx, atom2_idx);
                lj_energy += calculateLJ(r_ij);
            }
        }
    }
    return lj_energy;
}

const arma::mat& LennardJones::getPositions() const {
    return positions;
}

const std::vector<int>& LennardJones::getAtomicNumbers() const {
    return A;
}