// Solving 1D Steady State Convection Diffusion Equation using C++
// - A script to set up and solve the 1D convection diffusion equation for conduction in a bar using UPWING SCHEME
// [Chapter 4]

#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <cmath>
#include "matplotlibcpp.h"

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;
namespace plt = matplotlibcpp;


// Function to print a vector
void printVector(const std::vector<double>& vec, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

// Function to print a matrix
void printMatrix(const std::vector<std::vector<double>>& mat, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : mat) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<double> solveLinearSystem(Matrix A, Vector b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Swap maximum row with current row
        std::swap(A[maxRow], A[i]);
        std::swap(b[maxRow], b[i]);

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Solve Ux = y using back substitution
    Vector x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}


void ConvectionDiffusion1D(double barLength, int nCells, double area, double cond, double Uvel, double rho, double cP, double tempLeft, double tempRight, double heatSourcePerVol, bool printSetup, bool printSolution, bool printOutput)
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Generating Mesh ..." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    // Calculating Coordinates of Cell Faces
    std::vector<double> xFaces(nCells+1);
    double step = barLength / nCells;
    for(int i = 0; i <= nCells; i++)
    {
        xFaces[i] = i * step;
    }

    // Calculate the coordinates of Cell Centroids
    std::vector<double> xCentroids(nCells);
    for(int i = 0; i < nCells; i++)
    {
        xCentroids[i] = 0.5 * (xFaces[i+1] + xFaces[i]);
    }

    // Calculate Length of each Cell
    std::vector<double> cellLength(nCells);
    for(int i = 0; i < nCells; i++)
    {
        cellLength[i] = xFaces[i+1] - xFaces[i];
    }

    // Calculate distance between Cell Centroids (Diffusion dT/dx)
    std::vector<double> dCentroid(nCells - 1);
    for(int i = 0; i < nCells-1; i++)
    {
        dCentroid[i] = xCentroids[i+1] - xCentroids[i];
    }
    //dCentroid[0] = 2 * (xCentroids[0] - xFaces[0]);
    //dCentroid[nCells] = 2 * (xFaces[nCells] - xCentroids[nCells-1]);

    // Left Boundary Condition distance is half of distance between cell Face and cell Centroid
    double dLeft = 2*(xCentroids[0] - xFaces[0]);

    // Right Boundary Condition distance is half of distance between cell Face and cell Centroid
    double dRight = 2*(xFaces.back() - xCentroids.back());

    // Append the Vector of Cell Centroids
    dCentroid.insert(dCentroid.begin(),dLeft);
    dCentroid.push_back(dRight);

    //Compute the Cell Volume
    std::vector<double> cellVolume(nCells);
    for(int i = 0; i < nCells; i++)
    {
        cellVolume[i] = cellLength[i] * area;
    }

    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Mesh Generated ..." << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    //=============================MATRIX Setup==============================================================//
        // Diffusive Flux per unit Area
    // NOTE: DA represents the diffusive flux per unit area at the cell faces 
    //

    std::vector<double> DA(nCells+1);
    for(int i = 0; i <= nCells; i++)
    {
        DA[i] = ( area * cond )/ dCentroid[i];
    }

    //Convective Flux
    //F = rho * Cp * U * A
    // rho and U are Vectors
    std::vector<double> flowVelocity(nCells+1,Uvel);
    std::vector<double> flowDensity(nCells+1,rho);

    std::vector<double> F(nCells+1, 0.0);
    for(int i = 0; i < nCells + 1; i++)
    {
        F[i] = flowDensity[i] * flowVelocity[i] * area * cP; 
    }

    //Calculate the Peclet Number 
    std::vector<double> Pe(nCells+1,0.0);
    for(int i = 0; i < nCells + 1; i++)
    {
        Pe[i] = F[i] / DA[i]; 
    }   


    // Calculate the Source Term
    std::vector<double> Sp(nCells,0.0);

    //Assign Sp values to Left and Right Boundary Cell
    Sp[0] = -(2*DA[0] + std::max(F[0],static_cast<double>(0)));
    Sp[nCells - 1] = -(2*DA[nCells] + std::max(-1 * F[nCells], static_cast<double>(0)));

    //Calculate Source Term
    std::vector<double> Su(nCells);
    for(int i = 0; i < nCells; i++)
    {
        Su[i] = heatSourcePerVol * cellVolume[i];
    }

    //Assign Heat Source to the Left and Right Boundary Condition
    Su[0] = Su[0] + tempLeft*(2*DA[0] + std::max(F[0],static_cast<double>(0)));
    Su[nCells-1] = Su[nCells-1] + tempRight*(2*DA[nCells] + std::max(-1* F[nCells],static_cast<double>(0)));

    //aL Coefficients
    std::vector<double> aL(nCells);
    for(int i=0; i < nCells; i++)
    {
        aL[i] = DA[i] + std::max(F[i],static_cast<double>(0)) ;
    }
    aL[0] = 0.0;

    //aR Coefficients
    std::vector<double> aR(nCells);
    for(int i=0; i < nCells; i++)
    {
        aR[i] = DA[i] + std::max(-F[i], static_cast<double>(0));
    }
    aR[nCells-1] = 0.0;

    //aP Coefficient
    std::vector<double> aP(nCells);
    for(int i = 0; i < nCells; i++)
    {
        aP[i] = aL[i] + aR[i] + F[i+1] - F[i] - Sp[i];
    }

    //===============================Assembling Matrix=========================================//
    std::vector<std::vector<double>> AMatrix(nCells,std::vector<double> (nCells,0.0));
    std::vector<double> BMatrix = Su;

    //Allocating the Coefficients accordingly
    for(int i = 0; i < nCells; i++)
    {
        AMatrix[i][i] = aP[i];
        if(i > 0)
        {
            AMatrix[i][i-1] = -aL[i];
        }
        if(i < nCells - 1)
        {
            AMatrix[i][i+1] = -aR[i];
        }

        BMatrix[i] = Su[i];
    }
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << " Matrix Setup Complete ..." << std::endl;
        std::cout << "------------------------------------------------" << std::endl;

    //=============================Print the Setup===============================================//
    if(printSetup)
    {
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << " Summary: Set Up" << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        
        //printVector(xCentroids, "xC");
        //printVector(dCentroid, "dC");
        printVector(DA, "DA");
        printVector(aL, "aL");
        printVector(aR, "aR");
        printVector(aP, "aP");
        printVector(Sp, "Sp");
        printVector(Su, "Su");
        printVector(Pe,"Peclet Number");
        printMatrix(AMatrix, "A Matrix");
        printVector(BMatrix, "B Matrix");        
    }

    //===========================Solve the Matrices================================================//

        std::cout << "------------------------------------------------" << std::endl;
        std::cout << " Solving ... " << std::endl;
        std::cout << "------------------------------------------------" << std::endl;   

    Vector TVector = solveLinearSystem(AMatrix,BMatrix);


    //===============================Printing the Result===========================================//
    if(printSolution)
    {
        printVector(TVector,"Result T: ");
    }

    if (printOutput)
    {
        std::cout << "------------------------------------------------" << std::endl;
        std::cout << " Plotting ..." << std::endl;
        std::cout << "------------------------------------------------" << std::endl;

        // Assuming these variables are defined and calculated earlier in your code
        //std::vector<double> xFaces, xCentroids, Tvector;
        //double tempLeft, tempRight, barLength, heatSourcePerVol, cond;

        // Append the boundary temperature values to the vector
        std::vector<double> xPlotting = {xFaces.front()};
        xPlotting.insert(xPlotting.end(), xCentroids.begin(), xCentroids.end());
        xPlotting.push_back(xFaces.back());

        std::vector<double> temperaturePlotting = {tempLeft};
        temperaturePlotting.insert(temperaturePlotting.end(), TVector.begin(), TVector.end());
        temperaturePlotting.push_back(tempRight);

        // NO ANALYTICAL SOLUTION 
        /*
        // Assemble the analytical solution for comparison
        std::vector<double> xAnalytical(100);
        std::vector<double> temperatureAnalytical(100);
        for (int i = 0; i < 100; ++i) {
            xAnalytical[i] = i * barLength / 99.0;
            temperatureAnalytical[i] = tempLeft + ((tempRight - tempLeft) * (xAnalytical[i]/barLength)) +
                                       (heatSourcePerVol/(2.0*cond)) * xAnalytical[i] *
                                       (barLength - xAnalytical[i]);
        }
        */

        // Configure the plot
        //plt::figure_size(620, 420);
        plt::named_plot("CFD", xPlotting, temperaturePlotting, "-o");
        //plt::named_plot("Analytical", xAnalytical, temperatureAnalytical, "--");
        
        plt::xlabel("x [m]");
        plt::ylabel("T [°C]");
        plt::xlim(xFaces.front(), xFaces.back());
        plt::legend();

        // Set font sizes and line widths
        plt::title("Temperature Distribution",{{"fontsize", "14"}});

        //plt::save("temperature_plot.png");
        plt::show();
        

    }

}


int main(int argc,char* argv[])
{
    double barLength = 5;
    int nCells = 50;  // Increased for better resolution
    double area = 0.2;
    double cond = 200;
    double UVel = -0.1;
    double cP = 1000;
    double rho = 1.0;
    double heatSourcePerVol = 1000.0;  // Increased for more noticeable effect
    double tempLeft = 100;
    double tempRight = 400;
    bool printSetup = false;  // Set to false to reduce console output
    bool printSolution = true;
    bool printOutput = true;

    ConvectionDiffusion1D(barLength, nCells, area, cond, UVel, rho, cP, tempLeft, tempRight, heatSourcePerVol, printSetup, printSolution, printOutput);

    return 0;
}
