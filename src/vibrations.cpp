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

// Constants for conversion
const double conversion = 5140.484532;

void analyze_vibrations(const arma::mat & hessian) 
{
    int n_atoms = atomic_masses.size();
    int dof = 3 * n_atoms;

    // Build mass vector
    arma::vec masses(dof);
    for (int i = 0; i < n_atoms; i++) 
    {
        masses.subvec(3 * i, 3 * i + 2).fill(atomic_masses[i]);
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
