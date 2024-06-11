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
    FTL_Lesson06

Description
    This solver will introduce us to custom classes and as well how we can use derived 
    classes in our solver
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ftlClass.H"
#include "derivedClass.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //Create a custom class instance of type ftlclass, called custom FTL
    ftlClass customFTL;

    //Return the inline function get() to give us the value of myInt
    Info << "Default value for customFTL is: " << customFTL.get() << endl;

    //Set this to a new value using the inline function set()
    customFTL.set(10);

    //Rrint the new value of customFTL stored in myInt
    Info << "New value for CustomFTL is: " << customFTL.get() <<endl; 

    label ftlValue = customFTL.doubleStoredValue();
    Info << "New value for CustomFTL is: " << ftlValue <<endl; 

    customFTL.meshCountFunction(mesh);
    Info << " The value now stored in myInt_ for customFTL is " <<customFTL.get() <<endl;

    //*****************************************************************//
    ftlDict myTransportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    dimensionedScalar nu 
    (
        "nu",
        dimViscosity,
        myTransportProperties.lookup("nu")
    ); 
    Info << "Created a viscosity scalar: " <<nu <<endl;

    myTransportProperties.printTokensInTheDict();

    //*****************************************************************//
    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
