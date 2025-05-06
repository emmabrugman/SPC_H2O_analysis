#include "spc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

// Constants in atomic units
const double SPC::SIGMA_O = 5.9819129309; // Bohr
const double SPC::EPSILON_O = 0.0002476859286; // Hartree
const double SPC::COULOMB_CONSTANT = 1.0; // Dimensionless in atomic units

const double SPC::KB = 0.474; // Hartree/Bohr^2
const double SPC::KA = 0.121; // Hartree/rad^2
const double SPC::R_OH_EQ = 1.912; // Bohr
const double SPC::THETA_HOH_EQ = 1.976; // radians (113.24 degrees)

const double SPC::CUTOFF = 18.8973; // Bohr

// Constructor: initialize positions and atomic numbers
SPC::SPC(const arma::mat& pos, const std::vector<int>& a)
    : positions(pos), A(a) {
    charges.resize(A.size());
    for (size_t i = 0; i < A.size(); ++i) {
        charges[i] = getPartialCharge(A[i]);
    }
}

// Set positions
void SPC::setPositions(const arma::mat& newPositions) {
    positions = newPositions;
}

// SPC/E partial charge based on atomic number
double SPC::getPartialCharge(int atomic_number) {
    if (atomic_number == 8) return -0.8476; // Oxygen
    if (atomic_number == 1) return 0.4238;  // Hydrogen
    throw std::invalid_argument("Unsupported atomic number for SPC/E model");
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
    double r6_term = std::pow(SIGMA_O / r_ij, 6);
    double r12_term = r6_term * r6_term;  // (σ/r)^12
    return 4 * EPSILON_O * (r12_term - r6_term);  // Standard 12-6 LJ formula
}

double SPC::getLJEnergy() const {
    double lj_energy = 0.0;
    for (size_t atom1_idx = 0; atom1_idx < A.size(); ++atom1_idx) {
        for (size_t atom2_idx = atom1_idx + 1; atom2_idx < A.size(); ++atom2_idx) {
            if (A[atom1_idx] == 8 && A[atom2_idx] == 8) {
                double r_ij = calculateDistance(atom1_idx, atom2_idx);
                if (r_ij < CUTOFF) {
                    double lj_contrib = calculateLJ(r_ij);
                    lj_energy += lj_contrib;
                }
            }
        }
    }
    return lj_energy;
}

// ====== Coulomb Energy ======
double SPC::getCoulombEnergy() const 
{
    double energy = 0.0;
    for (size_t i = 0; i < A.size(); ++i) 
    {
        for (size_t j = i + 1; j < A.size(); ++j) 
        {
            double q1 = charges[i];
            double q2 = charges[j];
            double r = calculateDistance(i, j);
            if (r > 0 && r < CUTOFF) {
                double contrib = q1 * q2 / r;
                energy += contrib;
            }
        }
    }
    return energy;
}

// ====== Angle Energy ======
double SPC::getAngleEnergy() const {
    double angle_energy = 0.0;

    // Number of water molecules
    int num_molecules = A.size() / 3;
    
    // Process each molecule
    for (int mol = 0; mol < num_molecules; ++mol) {
        // Get atom indices for this molecule
        int oxygen_idx = mol * 3;       // First atom in the molecule
        int h1_idx = mol * 3 + 1;       // Second atom in the molecule
        int h2_idx = mol * 3 + 2;       // Third atom in the molecule
        
        // Get positions of oxygen and two hydrogens
        arma::vec pos_O = positions.col(oxygen_idx);
        arma::vec pos_H1 = positions.col(h1_idx);
        arma::vec pos_H2 = positions.col(h2_idx);

        // Compute vectors r_OH1 and r_OH2
        arma::vec r_OH1 = pos_H1 - pos_O;
        arma::vec r_OH2 = pos_H2 - pos_O;

        // Compute magnitudes
        double mag_OH1 = arma::norm(r_OH1, 2);
        double mag_OH2 = arma::norm(r_OH2, 2);

        // Compute cosine of the angle using dot product
        double cos_theta = arma::dot(r_OH1, r_OH2) / (mag_OH1 * mag_OH2);

        // Clamp cos_theta to [-1, 1] to avoid numerical errors
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

        // Compute the angle in radians
        double theta = std::acos(cos_theta);

        // Compute harmonic angle energy for this molecule
        double delta_theta = theta - THETA_HOH_EQ;
        angle_energy += 0.5 * KA * delta_theta * delta_theta;
    }

    return angle_energy;
}

// ====== Bond Energy ======
double SPC::getBondEnergy() const {
    double bond_energy = 0.0;

    // Check if the number of atoms is divisible by 3
    if (A.size() % 3 != 0) {
        throw std::runtime_error("Expected the number of atoms to be divisible by 3 (O, H, H pattern).");
    }

    // Number of water molecules
    int num_molecules = A.size() / 3;
    
    // Process each molecule
    for (int mol = 0; mol < num_molecules; ++mol) {
        // Get atom indices for this molecule
        int oxygen_idx = mol * 3;       // First atom in the molecule
        int h1_idx = mol * 3 + 1;       // Second atom in the molecule
        int h2_idx = mol * 3 + 2;       // Third atom in the molecule
        
        // Verify that we have O-H-H pattern
        if (A[oxygen_idx] != 8 || A[h1_idx] != 1 || A[h2_idx] != 1) {
            throw std::runtime_error("Expected O-H-H pattern for each water molecule.");
        }
        
        // Calculate bond energy for O-H1 bond
        double l1 = calculateDistance(oxygen_idx, h1_idx);
        double delta1 = l1 - R_OH_EQ;
        bond_energy += 0.5 * KB * delta1 * delta1;
        
        // Calculate bond energy for O-H2 bond
        double l2 = calculateDistance(oxygen_idx, h2_idx);
        double delta2 = l2 - R_OH_EQ;
        bond_energy += 0.5 * KB * delta2 * delta2;
    }

    return bond_energy;
}

// ====== Total Energy =====
double SPC::getTotalEnergy() const {
    return getLJEnergy() + getCoulombEnergy() + getBondEnergy() + getAngleEnergy();
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