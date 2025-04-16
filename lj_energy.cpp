#include "lj_energy.h"
#include <iostream>
#include <armadillo>
#include <cmath>

// Constructor
LennardJones::LennardJones(double eps, double sig) : epsilon(eps), sigma(sig) {}

// Calculate the total energy
double LennardJones::calculate_total_energy(const arma::mat & coords) 
{
    double total_energy = 0.0;

    // Check if the coords matrix is empty or has an incorrect number of columns
    if (coords.is_empty() || coords.n_cols != 3) 
    {
        std::cerr << "Error: the coordinates matrix is empty / has an invalid shape" << std::endl;
        exit(1);
    }

    // Loop through atom pairs
    for (size_t i = 0; i < coords.n_rows; i++) 
    {
        for (size_t j = i + 1; j < coords.n_rows; j++) 
        {
            // Calculate the vector between atoms i and j
            double r_ij = arma::norm(coords.row(i) - coords.row(j));

            // Check for division by zero error
            if (r_ij == 0.0) 
            {
                std::cerr << "Error: rij is zero" << std::endl;
                exit(1);
            }

            // Lennard Jones potential formula
            total_energy += epsilon * (pow(sigma / r_ij, 12) - 2 * pow(sigma / r_ij, 6));
        }
    }
    return total_energy;
}


