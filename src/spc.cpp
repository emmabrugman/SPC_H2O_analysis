#include "spc.h"
#include <stdexcept>
#include <cmath>

// Constants
const double SPC::SIGMA_O = 3.166; // Å
const double SPC::EPSILON_O = 0.650; // kJ/mol
const double SPC::COULOMB_CONSTANT = 332.0637; // kJ·Å/(mol·e²)

// Constructor: initialize positions and atomic numbers
SPC::SPC(const arma::mat& pos, const std::vector<int>& a)
    : positions(pos), A(a) {
    charges.resize(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        charges[i] = getPartialCharge(A[i]);
    }
}

// SPC partial charge based on atomic number
double SPC::getPartialCharge(int atomic_number) {
    if (atomic_number == 6) return -0.82; // Oxygen
    if (atomic_number == 1) return 0.41;  // Hydrogen
    throw std::runtime_error("Unknown atomic number: " + std::to_string(atomic_number));
}

// Compute Euclidean distance between two atoms
double SPC::calculateDistance(size_t first_atom_index, size_t second_atom_index) const {
    arma::vec diff = positions.col(first_atom_index) - positions.col(second_atom_index);
    return arma::norm(diff, 2);
}

// ====== Lennard-Jones Potential ======
double SPC::calculateR6Term(double r) const {
    return std::pow(SIGMA_O / r, 6);
}

double SPC::calculateR12Term(double r6_term) const {
    return r6_term * r6_term;
}

double SPC::calculateLJ(double r_ij) const {
    double r6_term = calculateR6Term(r_ij);
    double r12_term = calculateR12Term(r6_term);
    return EPSILON_O * (r12_term - 2 * r6_term);
}

double SPC::getLJEnergy() const {
    double lj_energy = 0.0;
    for (size_t atom1_idx = 0; atom1_idx < A.size(); ++atom1_idx) {
        for (size_t atom2_idx = atom1_idx + 1; atom2_idx < A.size(); ++atom2_idx) {
            // Apply LJ only if both atoms are oxygen (A=6)
            if (A[atom1_idx] == 6 && A[atom2_idx] == 6) {
                double r_ij = calculateDistance(atom1_idx, atom2_idx);
                lj_energy += calculateLJ(r_ij);
            }
        }
    }
    return lj_energy;
} 

// ====== Coulomb Potential ======
double SPC::getCoulombEnergy() const {
    double energy = 0.0;
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = i + 1; j < A.size(); ++j) {
            double q1 = charges[i]; // Use pre-computed charges
            double q2 = charges[j];
            double r = calculateDistance(i, j); // Use SPC's distance method
            if (r > 0) { // Prevent division by zero
                energy += COULOMB_CONSTANT * q1 * q2 / r;
            }
        }
    }
    return energy;
}

// ====== Total Energy =====
double SPC::getTotalEnergy() const {
    return getLJEnergy() + getCoulombEnergy();
}

// ====== Getters =====
const arma::mat& SPC::getPositions() const {
    return positions;
}

const std::vector<int>& SPC::getAtomicNumbers() const {
    return A;
}

const std::vector<double>& SPC::getCharges() const {
    return charges;
}