#pragma once

#include <armadillo>
#include <vector>
#include <cmath>

class SPC {
private:
    arma::mat positions; // Matrix of atomic positions
    std::vector<int> A; // Vector of atomic numbers
    std::vector<double> charges; // Vector of partial charges

public:
    // Constants
    static const double SIGMA_O; 
    static const double EPSILON_O; 
    static const double COULOMB_CONSTANT; 
    static const double CUTOFF; 
    static const double KB; // Bond stretching constant (Hartree/Bohr^2)
    static const double KA; // Angle bending constant (Hartree/rad^2)
    static const double R_OH_EQ; // Equilibrium O-H bond length (Bohr)
    static const double THETA_HOH_EQ; // Equilibrium H-O-H angle (radians)

    // Constructor: initialize positions, atomic numbers, and charges
    SPC(const arma::mat& pos, const std::vector<int>& a);

    // SPC/E partial charge based on atomic number
    double getPartialCharge(int atomic_number);

    // Compute the total LJ energy
    double getLJEnergy() const;

    // Compute the total Coulomb energy
    double getCoulombEnergy() const;

    // Compute the polarization correction energy (SPC/E specific)
    double getPolarizationEnergy() const;

    // Compute the total energy (LJ + Coulomb + Polarization)
    double getTotalEnergy() const;

    // Getters
    const arma::mat& getPositions() const;
    const std::vector<int>& getAtomicNumbers() const;
    const std::vector<double>& getCharges() const;

    // Setter for positions (for force calculations)
    void setPositions(const arma::mat& newPositions);

    // Calculates (sigma/r)^6
    double calculateR6Term(double r) const;

    // Calculates ((sigma/r)^6)^2
    double calculateR12Term(double r6_term) const;

    // Calculates Euclidean distance between two atoms
    double calculateDistance(size_t first_atom_index, size_t second_atom_index) const;

    // Calculates LJ potential for a given distance
    double calculateLJ(double r_ij) const;
};