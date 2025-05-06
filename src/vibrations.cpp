#include <armadillo>
#include <vector>
#include <cmath>
#include <iostream>
#include "vibrations.h"

// Atomic masses (amu)
static const std::vector<double> atomic_masses = 
{
    15.999, // O
    1.008,  // H
    1.008   // H
};

#include <armadillo>
#include <vector>
#include <cmath>
#include <iostream>
#include "vibrations.h"

// Constants for conversion
const double conversion = 5140.484532;

void analyze_vibrations(const arma::mat & hessian) 
{
    // Calculate dimensions
    int dof = hessian.n_rows;
    int n_atoms = dof / 3;
    
    // For multiple water molecules, we assume they are arranged as O, H, H, O, H, H, etc.
    int n_molecules = n_atoms / 3;
    
    std::cout << "Analyzing vibrational modes for " << n_molecules << " water molecule(s) "
              << "(" << n_atoms << " atoms, " << dof << " degrees of freedom)" << std::endl;

    // Build mass vector based on known water molecule pattern (O, H, H, O, H, H, ...)
    arma::vec masses(dof);
    for (int i = 0; i < n_atoms; i++) {
        if (i % 3 == 0) { // Oxygen (first atom in each molecule)
            masses.subvec(3*i, 3*i+2).fill(15.999);
        } else { // Hydrogen (second and third atoms in each molecule)
            masses.subvec(3*i, 3*i+2).fill(1.008);
        }
    }

    // Mass-weight the Hessian
    arma::mat mw_hessian(dof, dof);
    for (int i = 0; i < dof; i++) 
    {
        for (int j = 0; j < dof; j++) 
        {
            mw_hessian(i, j) = hessian(i, j) / std::sqrt(masses(i) * masses(j));
        }
    }

    // Diagonalize
    arma::vec eigvals;
    arma::mat eigvecs;
    arma::eig_sym(eigvals, eigvecs, mw_hessian);

    std::cout << "Eigenvalues (w^2 in a.u.):\n";
    eigvals.t().print();

    std::cout << "\nFrequencies (cm^-1):\n";
    for (double eig : eigvals) 
    {
        if (eig < 0) 
        {
            std::cout << "  Imaginary (" << eig << ")\n";
        } else {
            double freq = std::sqrt(eig) * conversion;
            std::cout << "  " << freq << " cm^-1\n";
        }
    }
}