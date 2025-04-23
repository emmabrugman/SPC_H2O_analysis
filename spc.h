#pragma once

#include <armadillo>
#include <vector>
#include <cmath>

class LennardJones {
private:
    arma::mat positions; // Matrix of atomic positions
    std::vector<int> A; // Vector of atomic numbers

public:
    // LJ Constants
    static const double SIGMA_O; // 3.166 Å
    static const double EPSILON_O; // 0.650 kJ/mol

    // Constructor: initialize positions and atomic numbers
    LennardJones(const arma::mat& pos, const std::vector<int>& a);

    // Compute energy
    double getLJEnergy() const;
    double getCoulombEnergy() const;
    double getTotalEnergy() const;

    // Getters
    const arma::mat& getPositions() const;
    const std::vector<int>& getAtomicNumbers() const;

    // Calculates (sigma/r)^6
    double calculateR6Term(double r) const;

    // Calculates ((sigma/r)^6)^2
    double calculateR12Term(double r6_term) const;

    // Calculates Euclidean distance between two atoms
    double calculateDistance(int first_atom_index, int second_atom_index) const;

    // Calculates LJ potential for a given distance
    double calculateLJ(double r_ij) const;
};
