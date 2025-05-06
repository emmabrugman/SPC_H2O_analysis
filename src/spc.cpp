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

// ====== Intra Energy ======
// Angle Energy +  Bond Energy

double SPC::getBondEnergy() const 
{
    double bond_energy = 0.0;

    // 1 oxygen (Z=8), 2 hydrogens (Z=1)
    int oxygen_idx = -1;
    std::vector<int> hydrogen_indices;
    for (size_t i = 0; i < A.size(); i++) 
    {
        if (A[i] == 8) 
        {
            oxygen_idx = i;
        } else if (A[i] == 1) {
            hydrogen_indices.push_back(i);
        }
    }

    // Check
    if (oxygen_idx == -1 || hydrogen_indices.size() != 2) 
    {
        throw std::runtime_error("Expected one O (Z=8) and two H (Z=1) atoms.");
    }

    // Harmonic bond energy
    for (int h_idx : hydrogen_indices) 
    {
        double l = calculateDistance(oxygen_idx, h_idx);  // bond length
        double delta = l - R_OH_EQ;                       // deviation from equilibrium
        bond_energy += 0.5 * KB * delta * delta;
    }

    return bond_energy;
}


// ====== Inter Energy ======
// LJ and Coulomb energy


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