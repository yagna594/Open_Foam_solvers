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
    ftlLesson09

Description
    Exploring the transport equation through defining a passive scalar transport

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    //solve the steady scalar transport equation using the solver specified in system/fvSolution dict 
    //Discretization of the individual term is specified in the system/fvSchemes
    solve
    (
        //Convective term advection of our passive scalar, T due to the velocity field
        //duT/dx
        fvm::div(phi,T)  
        
        //fvc::div(u) == fvc::div(phi)
        //[M] [O] = [B]
        //fvm:: generate coefficients for [M]
        // fvc:: goes in [B]

        //Diffiusive term - diffusion of our passive scalar, T, due to its own gradient and 
            //Proportionality constants DT, i.e. diffusion coefficient
        
        // DT * d2T/dx2 
        - fvm::laplacian(DT,T) //[m2/s]*[1/m2]*[kg/m3] == [kg m3/s] <=> flux of T

        //Source term
        //** fvOptions(T)
    );

    volScalarField result
    (
        IOobject
        (
            "result",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T
    );
    result.write();
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
