#include "build_system.h"
#include "spc.h"
#include "forces.h"
#include "hessian.h"
#include "vibrations.h"
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

    // Display number of water molecules
    int num_molecules = atomic_system.getNumAtoms() / 3;
    std::cout << "\n======== System Information ========" << std::endl;
    std::cout << "Number of atoms: " << atomic_system.getNumAtoms() << std::endl;
    std::cout << "Number of water molecules: " << num_molecules << std::endl;

    // Compute and display energies
    double lj_energy = spc_calculator.getLJEnergy();
    double coulomb_energy = spc_calculator.getCoulombEnergy();
    double bond_energy = spc_calculator.getBondEnergy();
    double angle_energy = spc_calculator.getAngleEnergy();
    double total_energy = spc_calculator.getTotalEnergy();

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\n======== Energy Calculations ========" << std::endl;
    std::cout << "Lennard-Jones Energy: " << lj_energy << " Hartree" << std::endl;
    std::cout << "Coulomb Energy: " << coulomb_energy << " Hartree" << std::endl;
    std::cout << "Bond Energy: " << bond_energy << " Hartree" << std::endl;
    std::cout << "Angle Energy: " << angle_energy << " Hartree" << std::endl;
    std::cout << "Total Energy: " << total_energy << " Hartree" << std::endl;
    
    std::cout << "\n======== Force Calculation ========" << std::endl;
    CalculateForce force_calculator(spc_calculator, atomic_system.getNumAtoms());
    arma::mat central_forces;
    double h_force = 0.0188973; // 0.01 Å in Bohr
    force_calculator.calculateCentralDifferenceForces(h_force, central_forces);
    
    // Output forces
    std::cout << "\nCentral Difference Forces (Hartree/Bohr):\n";
    std::cout << central_forces << std::endl;
    
    std::cout << "\n======== Hessian Calculation ========" << std::endl;
    CalculateHessian hessian_calculator(spc_calculator, force_calculator, atomic_system.getNumAtoms());
    arma::mat hessian;
    double h_hessian = 0.0188973; // Same step size as forces for consistency
    hessian_calculator.calculateCentralDifferenceHessian(h_hessian, hessian);
    
    // Output Hessian
    std::cout << "\nNumerical Hessian (Hartree/Bohr^2):\n";
    std::cout << hessian << std::endl;

    // Analyze vibrations
    std::cout << "\n======== Vibrational Analysis ========" << std::endl;
    analyze_vibrations(hessian);
    
    return 0;
}