#include "build_system.h"
#include "spc.h"
#include "forces.h"
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

    // Initialize SPC calculator
    SPC spc_calculator(atomic_system.getPositions(), atomic_system.getAtomicNumbers());

    // Compute and display energies
    double lj_energy = spc_calculator.getLJEnergy();
    double coulomb_energy = spc_calculator.getCoulombEnergy();
    double total_energy = spc_calculator.getTotalEnergy();

    std::cout << std::fixed << std::setprecision(5);
    std::cout << " \n======== Initial Energy Calculations -======== \n" << std::endl;
    std::cout << "Lennard-Jones Energy: " << lj_energy << " kJ/mol" << std::endl;
    std::cout << "Coulomb Energy: " << coulomb_energy << " kJ/mol" << std::endl;
    std::cout << "Total Energy (LJ + Coulomb): " << total_energy << " kJ/mol" << std::endl;

    CalculateForce force_calculator(spc_calculator, atomic_system.getNumAtoms());
    arma::mat central_forces;
    double h = 1e-5; // Small displacement for numerical derivative
    force_calculator.calculateCentralDifferenceForces(h, central_forces);
    
    // Output forces
    std::cout << "\nCentral Difference Forces (kJ/(mol·Å)):\n";
    std::cout << central_forces << std::endl;

return 0;
}