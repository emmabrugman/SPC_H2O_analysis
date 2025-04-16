#pragma once

#include "atom.h"
#include <armadillo>
#include <vector>

class BuildSystem  {
private:
    arma::mat positions; // Matrix of size 3 x N: each column holds {x, y, z}.
    std::vector<int> A;  // Vector storing the atomic numbers.
    int num_atoms;       // Total number of atoms in the system.

public:
    // Reads atomic data from file
    void readFromFile(const std::string& filename); 

    // Adds atom to the system (matrix)
    void addAtom(const Atom& atom); 

    // Returns a const reference to the positions
    const arma::mat& getPositions() const; 

    // Returns a const reference to the atomic numbers
    const std::vector<int>& getAtomicNumbers() const;

    // Returns the number of atoms in the system
    int getNumAtoms() const { return num_atoms; } 
};