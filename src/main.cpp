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

    // Initialize SPC/E calculator
    SPC spc_calculator(atomic_system.getPositions(), atomic_system.getAtomicNumbers());

    // Compute and display energies
    double lj_energy = spc_calculator.getLJEnergy();
    double coulomb_energy = spc_calculator.getCoulombEnergy();
    double polarization_energy = spc_calculator.getPolarizationEnergy();
    double total_energy = spc_calculator.getTotalEnergy();

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\n======== Initial Energy Calculations ========" << std::endl;
    std::cout << "Lennard-Jones Energy: " << lj_energy << " Hartree" << std::endl;
    std::cout << "Coulomb Energy: " << coulomb_energy << " Hartree" << std::endl;
    std::cout << "Polarization Energy: " << polarization_energy << " Hartree" << std::endl;
    std::cout << "Total Energy (LJ + Coulomb + Polarization): " << total_energy << " Hartree" << std::endl;
    
    std::cout << "\n======== Force Calculation ========" << std::endl;
    CalculateForce force_calculator(spc_calculator, atomic_system.getNumAtoms());
    arma::mat central_forces, analytical_forces;
    double h_force = 1.88972612457e-8; // 1e-8 Å converted to Bohr
    force_calculator.calculateCentralDifferenceForces(h_force, central_forces);
    
    // Output forces
    std::cout << "\nCentral Difference Forces (Hartree/Bohr):\n";
    std::cout << central_forces << std::endl;
    
    return 0;
}