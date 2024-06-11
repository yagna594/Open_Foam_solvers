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
    FTL-Lesson02

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Initialise of case
    #include "setRootCase.H"
    // This creates the time system for is (instance called runTime)
    #include "createTime.H"
    // This creates the time system for is (instance called runTime)
    #include "createMesh.H"

    dictionary customDict;
    const word dictName("flowthermolabProperties");

    IOobject dictIO
    {
        dictName, //This is the name of the file
        mesh. time() . constant (), //Path to the file
        mesh, // Reference to the mesh needed by the constructor
        IOobject :: MUST_READ // Means that the file must exist at runtime
    };

    if (!dictIO.typeHeaderOk<dictionary>(true))
        FatalErrorIn(args. executable()) << "Cannot open the custom FTL dictionary "
            << dictName << exit(FatalError);

    customDict = IOdictionary(dictIO);

    word someUsername;
    customDict.lookup("username")>> someUsername;


    scalar offsetValue(customDict.lookupOrDefault<scalar>("offsetValue",5.0));

    bool incompressibleFlag(customDict.lookupOrDefault<Switch>("incompressibleFlag",true));

    List< scalar> inputValues (customDict.lookup("inputValues"));

    HashTable<vector,word> inputHashTable(customDict.lookup("inputHashTable"));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<<nl<< "Looked up the following information:" <<nl <<nl
        <<"Solver run by username: " << someUsername <<nl
        <<"Your offset value is : "<<offsetValue <<nl
        << "You are running the simulation as incompressible: " << incompressibleFlag << nl
        << "With input values: " << inputValues << nl
        << "With custom hashTable: " << inputHashTable << endl;

    // Create the output path directory
    fileName outputDirectory = mesh. time().path()/"postProcessing";
    // Create the directory
    mkDir(outputDirectory);

    // File pointer to direct the output to
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr. reset (new OFstream(outputDirectory/"simulationProperties.dat"));

    // Write to the output file
    outputFilePtr() << "# Simulation Information" << endl;
    outputFilePtr() << "# Simulation was run by username:" << someUsername << endl;

    inputHashTable.insert("pressure field",vector(13,14,15));
    outputFilePtr() << inputHashTable <<endl;

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
