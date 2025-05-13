# Vibrational Analysis of Water using SPC Force Field Models

This repository contains the final project paper and slides, as well as the supporting code for vibrational analysis of water molecules using Single Point Charge (SPC) force field model variants, specifically SPC/E and SPC/Fw.

## Paper Abstract

This study presents a computational vibrational analysis of water molecules using Single Point Charge (SPC) force field model variants (SPC/E and SPC/Fw). We developed a C++ program to simulate the atomic system of a water molecule, compute intramolecular and intermolecular energies, and derive forces and the Hessian matrix via central difference approximations. Vibrational modes were obtained by diagonalizing the mass-weighted Hessian, yielding frequencies in wavenumber.

## Getting Started

1. Clone the repository:
   ```
   git clone git@github.com:emmabrugman/SPC_H2O_analysis.git
   cd SPC_H2O_analysis
   ```

2. Check out the appropriate branch:
   - For SPC/E model: `git checkout spc/e`
   - For SPC/Fw model: `git checkout spc/fw`

### SPC/FW Version

#### Repository Structure

```
├── build.sh              # Script to build the project
├── clean.sh              # Script to clean build files
├── CMakeLists.txt        # Main CMake configuration
├── input                 # Directory containing input files
    ├── h2o_HBOND5A.txt   # Sample water molecule geometry for 3 molecules 5 A apart
    ├── h2o_HBOND2-2A.txt # Sample water molecule geometry for 3 molecules 2.2 A apart
    ├── h2o_3.txt         # Sample water molecule geometry for 3 molecules 3 A apart
│   └── H2O.txt           # Sample water molecule geometry for single molecule
├── interactive.sh        # Script to launch Docker container
├── README.md             # This file
└── src                   # Source code directory
    ├── atom.h            # Atom class definition
    ├── build_system.cpp  # System building implementation
    ├── build_system.h    # System building header
    ├── CMakeLists.txt    # Source-specific CMake config
    ├── forces.cpp        # Force calculation implementation
    ├── forces.h          # Force calculation header
    ├── hessian.cpp       # Hessian matrix implementation
    ├── hessian.h         # Hessian matrix header
    ├── main.cpp          # Main program entry point
    ├── spc.cpp           # SPC model implementation
    ├── spc.h             # SPC model header
    ├── vibrations.cpp    # Vibrational analysis implementation
    └── vibrations.h      # Vibrational analysis header
└── documents             # files for submission
    ├── SPC_PPT.pdf       # Powerpoint used for final presentation
    └── SPC_Report.pdf    # Final report
```

### SPC/E Version

#### Repository Structure

```
├── build.sh              # Script to build the project
├── clean.sh              # Script to clean build files
├── CMakeLists.txt        # Main CMake configuration
├── input                 # Directory containing input files
│   └── H2O.txt           # Sample water molecule geometry
├── interactive.sh        # Script to launch Docker container
├── README.md             # This file
└── src                   # Source code directory
    ├── atom.h            # Atom class definition
    ├── build_system.cpp  # System building implementation
    ├── build_system.h    # System building header
    ├── CMakeLists.txt    # Source-specific CMake config
    ├── forces.cpp        # Force calculation implementation
    ├── forces.h          # Force calculation header
    ├── hessian.cpp       # Hessian matrix implementation
    ├── hessian.h         # Hessian matrix header
    ├── main.cpp          # Main program entry point
    ├── spc.cpp           # SPC model implementation
    ├── spc.h             # SPC model header
    ├── vibrations.cpp    # Vibrational analysis implementation
    └── vibrations.h      # Vibrational analysis header
```

### Running the Code

1. Start the Docker container environment:
   ```
   ./interactive.sh
   ```

2. Build the project:
   ```
   ./build.sh
   ```

3. Run the executable:
   ```
   ./build/src/spc
   ```

4. When prompted, enter the path to the input file:
   ```
   Enter the path to the atomic structure file: input/{input_file}.txt
   ```
