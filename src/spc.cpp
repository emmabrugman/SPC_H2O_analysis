// spc.cpp
#include "spc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

// Constants
const double SPC::SIGMA_O = 3.16555789; // Å
const double SPC::EPSILON_O = 0.155213; // kcal/mol
const double SPC::COULOMB_CONSTANT = 332.0637; // kcal·Å/(mol·e²)
const double SPC::POLARIZATION_ENERGY = 1.247; // kcal/mol (equivalent to 5.22 kJ/mol)
const double CUTOFF = 10.0; // Å

// Constructor: initialize positions and atomic numbers
SPC::SPC(const arma::mat& pos, const std::vector<int>& a)
    : positions(pos), A(a) {
    charges.resize(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        charges[i] = getPartialCharge(A[i]);
    }
}

// SPC/E partial charge based on atomic number
double SPC::getPartialCharge(int atomic_number) {
    if (atomic_number == 8) return -0.8476; // Oxygen
    if (atomic_number == 1) return 0.4238;  // Hydrogen
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
    int pair_count = 0;
    for (size_t atom1_idx = 0; atom1_idx < A.size(); ++atom1_idx) {
        for (size_t atom2_idx = atom1_idx + 1; atom2_idx < A.size(); ++atom2_idx) {
            if (A[atom1_idx] == 8 && A[atom2_idx] == 8) {
                double r_ij = calculateDistance(atom1_idx, atom2_idx);
                if (r_ij < CUTOFF) {
                    double lj_contrib = calculateLJ(r_ij);
                    lj_energy += lj_contrib;
                    pair_count++;
                }
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
            double q1 = charges[i];
            double q2 = charges[j];
            double r = calculateDistance(i, j);
            if (r > 0 && r < CUTOFF) {
                energy += COULOMB_CONSTANT * q1 * q2 / r;
            }
        }
    }
    return energy;
}

// ====== Polarization Correction (SPC/E specific) ======
double SPC::getPolarizationEnergy() const {
    // Count the number of water molecules
    // For SPC/E, each O atom corresponds to one water molecule
    int water_count = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] == 8) water_count++;
    }
    
    // The polarization energy is a constant per water molecule
    return water_count * POLARIZATION_ENERGY;
}

// ====== Total Energy =====
double SPC::getTotalEnergy() const {
    return getLJEnergy() + getCoulombEnergy() + getPolarizationEnergy();
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