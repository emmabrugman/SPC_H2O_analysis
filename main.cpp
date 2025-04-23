#include "build_system.h"
//#include "lennard_jones.h"
#include "spc.h"
#include "coulomb.h"
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char** argv) {
    std::string filename;

    if (argc > 1) {
        filename = argv[1];
    } else {
        std::cout << "Enter the path to the atomic structure file: ";
        std::cin >> filename;
    }

    // Create and populate the atomic system
    BuildSystem atomic_system;
    try {
        atomic_system.readFromFile(filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error reading file: " << e.what() << std::endl;
        return 1;
    }

    // Initialize Lennard-Jones calculator
    LennardJones lj_calculator(atomic_system.getPositions(), atomic_system.getAtomicNumbers());

    // Compute and display LJ energy
    double lj_energy = lj_calculator.getLJEnergy();

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Lennard-Jones Energy: " << lj_energy << " kJ/mol" << std::endl;
    std::cout << "Coulomb Energy: " << lj_calculator.getCoulombEnergy() << " kJ/mol" << std::endl;
    std::cout << "Total Energy: " << lj_calculator.getTotalEnergy() << " kJ/mol" << std::endl;

    return 0;
}