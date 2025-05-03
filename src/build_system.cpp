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
    const double angstrom_to_bohr = 1.88972612457;

    // Read each atom's data (atomic number, x, y, z)
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
        file >> A >> x >> y >> z;
        double x_bohr = x * angstrom_to_bohr;
        double y_bohr = y * angstrom_to_bohr;
        double z_bohr = z * angstrom_to_bohr;
        std::cout << "Atom " << atom_index << ": A=" << A << ", (x,y,z) in Bohr=(" << x_bohr << "," << y_bohr << "," << z_bohr << ")\n";
        addAtom(Atom(A, x_bohr, y_bohr, z_bohr));
    }
}

void BuildSystem::addAtom(const Atom& atom) {
    positions.insert_cols(positions.n_cols, atom.position);
    A.push_back(atom.A);
}

const arma::mat& BuildSystem::getPositions() const {
    return positions;
}

const std::vector<int>& BuildSystem::getAtomicNumbers() const {
    return A;
}