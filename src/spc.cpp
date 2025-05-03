#include "spc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

// Constants in atomic units
const double SPC::SIGMA_O = 5.982037444913; // Bohr
const double SPC::EPSILON_O = 0.000247347607; // Hartree
const double SPC::COULOMB_CONSTANT = 1.0; // Dimensionless in atomic units
const double SPC::POLARIZATION_ENERGY = -0.0019872; // Hartree
const double CUTOFF = 18.8973; // Bohr

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
    else {
        throw std::invalid_argument("Unsupported atomic number for SPC/E model");
    }
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

// ====== Coulomb Energy ======
double SPC::getCoulombEnergy() const {
    double energy = 0.0;
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = i + 1; j < A.size(); ++j) {
            double q1 = charges[i];
            double q2 = charges[j];
            double r = calculateDistance(i, j);
            if (r > 0 && r < CUTOFF) {
                double contrib = q1 * q2 / r;
                energy += contrib; // Added semicolon here
            }
        }
    }
    return energy;
}

// ====== Polarization Energy ======
double SPC::getPolarizationEnergy() const {
    size_t num_oxygen = 0;
    for (size_t i = 0; i < A.size(); ++i) {
        if (A[i] == 8) { // Oxygen
            num_oxygen++;
        }
    }
    return num_oxygen * POLARIZATION_ENERGY; // Per water molecule
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