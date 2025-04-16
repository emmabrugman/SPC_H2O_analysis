#pragma once

#include <armadillo>

struct Atom
{
    int A; // Atomic Z-number
    arma::vec3 position; // Armadillo vector of size 3 for {x, y, z}

    // Constructor: initializes an atom with atomic number and position
    Atom(int A, double x, double y, double z) 
        : A(A), position(arma::vec3({x, y, z})) {}
};