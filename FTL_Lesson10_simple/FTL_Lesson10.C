/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    FTL_Lesson10

Description
    This will be a solver to implement a light weight version of the SIMPLE
    algorithms, much like what is available in "simpleFoam". However this solver 
    demonstrates how the code works and the optimization or convergence check are
    ignored in this code.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    
    //Relaxation factor
    scalar alpha;
    fvSolution.lookup("alpha") >> alpha;
    //Cell containing reference pressure
    scalar pRefCell;
    fvSolution.lookup("pRefCell") >> pRefCell;
    //Reference pressure value
    scalar pRefValue;
    fvSolution.lookup("pRefValue") >> pRefValue;

    Info << "Read the following parameters:" <<endl;
    Info <<"Relaxation Factor alpha" << alpha <<endl;
    Info << " Index of cell for reference pressure: " <<pRefCell;
    Info << "With a refernce pressure value: "<< pRefValue <<endl;
    
    while (runTime.loop())
    {
        Info << nl<< "Iteration: " << runTime.timeName() <<endl;

        //Define the momentum equation:
        fvVectorMatrix UEqn
        (
            fvm::div(phi,U) - fvm::laplacian(nu,U) == -fvc::grad(p)
        );

        //Solve momentum equation for the current value of pressure
        UEqn.solve();

        volScalarField A = UEqn.A();
        volVectorField H = UEqn.H();
        //AU = H - nabla(P)
        //U = H/A - (1/A)* nabla(P)

        //nabla(U) = 0

        //nabla((1/A)*nabla(p))= nabla(H/A)
        volScalarField A_inv = 1.0/A; //Computing the inverse of A
        surfaceScalarField A_inv_flux = fvc :: interpolate(A_inv); //Interpolating inverse of A for laplacian operator

        volVectorField HbyA = A_inv * H; 

        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv_flux,p) == fvc::div(HbyA)
        );
        
        //Set the reference pressure for the equation
        pEqn.setReference(pRefCell,pRefValue);

        //Solve the pressure correction equation
        pEqn.solve();

        //Explicit under relaxation of the pressure equation
        p = alpha*p+(1.0-alpha)*p_old;

        //update the velocity field with the newly computed pressure field
        U = (A_inv*H) -(A_inv*fvc::grad(p));
        phi = fvc::interpolate(U) & mesh.Sf();

        //Update the boundary conditions for both p and U fields
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        p_old = p;

        runTime.timeName();
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
