#include "build_system.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

void BuildSystem::readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    file >> num_atoms;

    std::string buffer;
    std::getline(file, buffer);

    int A;
    double x, y, z;

    // Read each atom's data (atomic number, x, y, z)
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
        file >> A >> x >> y >> z; 

        addAtom(Atom(A, x, y, z));
    }
}

void BuildSystem::addAtom(const Atom& atom) {
    // Insert the atoms posisiton to the positions matrix
    positions.insert_cols(positions.n_cols, atom.position);
    // Insert the atomic number to the A vector
    A.push_back(atom.A);
}

const arma::mat& BuildSystem::getPositions() const {
    return positions;
}

const std::vector<int>& BuildSystem::getAtomicNumbers() const {
    return A;
}