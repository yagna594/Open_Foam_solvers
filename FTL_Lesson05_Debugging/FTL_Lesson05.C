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
    FTL-Lesson05

Description
    The solver shows us how we can start manipulating basic openfoam Velocity and 
    pressure fields. We will start to look at how oue fields are dimensioned and 
    read ito OF:
    We will then modify our domain ro apply a time varying pressure.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar calculatePressure(scalar t, vector x, vector x0, scalar scale);

// ************************************************************************* //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Reading in the transport propertoes" << nl <<endl;

    IOdictionary  transportProperties
    (
        IOobject
        (
            "transportProperties", //dictionary name
            runTime.constant(),    //where the dict file is stored
            mesh,
            IOobject::MUST_READ_IF_MODIFIED, //properties should be re read back during run time
            IOobject::NO_WRITE   // no need to write anything to transportproperties
        )
    );

    dimensionedScalar nu
    ("nu",dimViscosity,transportProperties.lookup("nu")); //dimensionSet(0,2,-1,0,0,0,0),

    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(), // time name is specified in controlDict like start time which starts at 0 folder
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE // we need to write as specified in controlDict as write time
        ),
        mesh //we need to initialize this field 
    );

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        ("U", runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),
        mesh
    );

    const vector originVector(0.05,0.05,0.05);
    const scalar rFarCell = max(mag(dimensionedVector("x0",dimLength,originVector)-mesh.C())).value();
    Info << nl<< "rFarCell is : " << rFarCell <<endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info << "\nStarting time loop\n" <<endl;

    while(runTime.loop())
    {
        Info <<"Time = " << runTime.timeName() <<nl <<endl;

        for (label cellI=0; cellI<mesh.C().size(); cellI++)
        {
            //calculate pressure is a own function
            p[cellI] = calculatePressure(runTime.time().value(),mesh.C()[cellI], originVector, rFarCell);
        }
        Info << nl<< "p[0] is : " << p[0] <<endl;
        U=fvc::grad(p)*dimensionedScalar("tmp",dimTime,1.0);
        Info << nl<< "U[0] is : " << U[0] <<endl;

        //units of Pressure = N/m2 ==kg/ms2
        //units of grad(pressure)= kg/m2s2
        //p* = p/rho = [kg/ms2]/[kg/m3] == [m2/s2]
        //units of grad(p*) == m/s2
        runTime.write();
    }

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //


//Define custonm pressure function
scalar calculatePressure(scalar t, vector x, vector x0, scalar scale)
{
    scalar r (mag(x-x0)/scale);

    scalar rR(1.0/(r+1e-12));

    scalar f(1.0);

    return Foam::cos(2.0*Foam::constant::mathematical::pi*f*t)*rR;
}

// ************************************************************************* //
