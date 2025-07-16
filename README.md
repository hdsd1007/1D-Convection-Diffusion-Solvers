# Project Overview

This repository provides C++ implementations designed to solve the one-dimensional (1D) steady-state convection-diffusion equation using finite volume methods. These solvers are fundamental in computational fluid dynamics (CFD) and heat transfer, offering numerical solutions for various physical phenomena, such as heat conduction in a bar or contaminant transport in a fluid.

The project includes three distinct solvers, each addressing a specific aspect or numerical scheme:

- `1D Diffusion Solver` (solve1DDiffusion.cpp): Focuses on pure diffusion problems, often representing heat conduction.

- `1D Convection-Diffusion Solver` (solve1DConvection_Diffusion.cpp): Combines both convection and diffusion phenomena, demonstrating their interplay.

- `1D Convection-Diffusion Upwind Solver` (solve1DConvection_DiffusionUpwind.cpp): Implements the Upwind Differencing Scheme to address numerical stability issues that can arise in convection-dominated flows, particularly important for preventing non-physical oscillations.

Mathematical Formulation
The generalized 1D steady-state convection-diffusion equation can be expressed as:

$$ \nabla \cdot (\rho c_p U \phi) = \nabla \cdot (\gamma \nabla \phi) + S $$

Where:

- $\phi$ is the dependent variable (e.g., temperature, concentration).

- $\rho$ is the fluid density.

- U is the fluid velocity.

- $\gamma$ is the diffusion coefficient (e.g., thermal conductivity for heat transfer, mass diffusivity for species transport).

- S is the source term.

**Upwind Differencing Scheme**: 
The solve1DConvection_DiffusionUpwind.cpp specifically utilizes the Upwind Differencing Scheme for the convective term. This scheme approximates the convective flux at a cell face using the value of 
phi from the upstream node (in the direction of flow). This approach inherently introduces numerical diffusion but is highly effective in maintaining solution stability and preventing non-physical oscillations, especially when the Peclet number ($$ Pe=
\frac{U * L}{\gamma} $$) is high, indicating convection-dominated flow.

The discrete form of the coefficients for the Upwind scheme is derived based on the direction of flow to ensure proper approximation of the convective term.

# Features

**Finite Volume Method**: All solvers are built upon the finite volume approach, discretizing the domain into a series of control volumes.

**Tri-diagonal Matrix Solver**: The resulting system of linear algebraic equations is solved efficiently using Gaussian Elimination with pivoting, implemented for tri-diagonal matrices.

**Flexible Boundary Conditions**: Supports Dirichlet boundary conditions for both left and right ends of the domain.

**Configurable Parameters**: Easily adjust physical properties such as bar length, number of cells, cross-sectional area, thermal conductivity, fluid velocity, density, specific heat capacity, and volumetric heat source within the source code.

**matplotlib-cpp Integration**: The solve1DConvection_DiffusionUpwind.cpp leverages the matplotlib-cpp library for direct plotting of temperature distributions, providing visual analysis of the numerical solutions.

# Getting Started

Prerequisites
To compile and run these solvers, you will need:

- A C++ compiler (e.g., g++, MinGW, MSVC).

- Python 3.x with matplotlib installed (for matplotlib-cpp visualization).

- The matplotlib-cpp library (available at https://github.com/lava/matplotlib-cpp). Follow its installation instructions, especially linking with Python, as this library is crucial for the plotting functionality in the Upwind solver.

# Compilation
Navigate to the project root directory in your terminal and compile the desired solver. For example, to compile the Upwind solver:

```bash
g++ "solve1DConvection_DiffusionUpwind.cpp" -o 1D_Convection_Diffusion_Upwind.exe -lmatplotlib-cpp -I/usr/include/python3.8 -lpython3.8
### (Adjust Python include path and library names as per your system setup and OS)
```
You can optionally create an executables folder to store your compiled binaries separately:

```bash
mkdir executables
g++ "solve1DConvection_DiffusionUpwind.cpp" -o executables/1D_Convection_Diffusion_Upwind.exe ...
```

Repeat this compilation process for solve1DDiffusion.cpp and solve1DConvection_Diffusion.cpp to create their respective executables.

Running the Solvers
Once compiled, you can run the executables from your terminal. If you placed them in an executables folder:

```bash
./executables/1D_Convection_Diffusion_Upwind.exe
The output, including the temperature distribution plot (for the Upwind solver), will be displayed.
```

Project Structure
.
├── solve1DDiffusion.cpp
├── solve1DConvection_Diffusion.cpp
├── solve1DConvection_DiffusionUpwind.cpp
└── README.md

# Other files in your local directory are ignored by .gitignore as per your setup.
Contributing
Feel free to fork this repository, open issues, or submit pull requests if you have suggestions for improvements or bug fixes.