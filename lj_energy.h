#ifndef LJ_ENERGY_H
#define LJ_ENERGY_H

#include <armadillo>
#include <cmath>

class LennardJones {
/*
Represents the Lennard-Jones Potential to help compute energy and force calculations 
*/
protected:
    double epsilon; //kcal/mol; Well depth or the strength of the interactions between atoms
    double sigma;   //angstroms; Distance between atoms where the interaction energy is at a minimum

public:
    // Initialize sigma and epsilon
    LennardJones(double eps, double sig);

    // Function to calculate the total energy of a system
    double calculate_total_energy(const arma::mat & coords);

};

#endif // LJ_ENERGY_H

