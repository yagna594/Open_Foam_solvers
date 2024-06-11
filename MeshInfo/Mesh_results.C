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
    MeshInfo
Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "The most recent time folder is "<<runTime.timeName()<<nl
        << "The mesh has " << mesh.C().size() << "cells."/*and "<< mesh.Cf().size()
        <<" internal faces. " */<< endl;
/*    
    for (label cellI =0; cellI <mesh.C().size();cellI++)
        if (cellI%50 == 0)
            Info << "Cell " << cellI << " this has a center at " << mesh.C()[cellI]<< endl;
    Info << endl;

    for (label faceI = 0; faceI <mesh.owner().size(); faceI++)
        if (faceI%100 ==0)
             Info << "Internal face "<<faceI << " with centre at " << mesh.Cf()[faceI]
             << " belongs to cell " << mesh.owner()[faceI]
             <<" has a neighbour " << mesh.neighbour()[faceI] << " "<< mesh.Cf()[mesh.neighbour()[faceI]] << endl;
    Info <<endl;

    forAll(mesh.boundaryMesh(),patchI)
        Info << " Patch "<< patchI << " = " <<mesh.boundaryMesh()[patchI].name() << " with " 
        << mesh.boundary()[patchI].Cf().size() << " faces. This starts at face " <<mesh.boundary()[patchI].start()<<endl;
    Info <<endl;
*/
    for (label patchI = 0; patchI <mesh.boundaryMesh().size(); patchI++)
    {
        label patchID(0);
        scalar sufce_Area = 0;  
        patchID =mesh.boundaryMesh().findPatchID(mesh.boundaryMesh()[patchI].name());
        forAll(mesh.boundary()[patchID].Sf(),faceI)
        {
            sufce_Area += mag(mesh.boundary()[patchID].Sf()[faceI]);
        }
        Info << "Surface Area of " << mesh.boundaryMesh()[patchI].name() << " = " <<sufce_Area <<endl;  
    }
    Info <<endl;

/*    label patchFaceI(0);
    forAll(mesh.boundaryMesh(),patchI)
        Info <<"Patch " << patchI << " has the face "<< patchFaceI << " adjacent to cell "
        << mesh.boundary()[patchI].patch().faceCells()[patchFaceI]
        << ". The normal vector of this is "<< mesh.boundary()[patchI].Sf()[patchFaceI]
        <<"and Surface area " <<mag(mesh.boundary()[patchI].Sf()[patchFaceI])<<endl;
    Info <<endl;

    const faceList fcs = mesh.faces();
    const List<point>& pts =mesh.points();
    const List<point>& centers = mesh.faceCentres();

    forAll(fcs, faceI)
        if(faceI%100==0)
        {
            if(faceI<mesh.Cf().size())
                Info <<"Internal face ";
            else 
            {
                forAll(mesh.boundary(), patchI)
                    if((mesh.boundary()[patchI].start()<= faceI) &&
                        (faceI < mesh.boundary()[patchI].start()+mesh.boundary()[patchI].Cf().size()))
                    {
                        Info << " Face is on patch " <<patchI << ", faceI ";
                        break;
                    }
            }
            Info << faceI << " With centre at " <<centers[faceI]
            << " have " <<fcs[faceI].size() << " Vertices: ";
            forAll(fcs[faceI], vertexI)
                Info <<" " <<pts[fcs[faceI][vertexI]];
            Info << endl; 
        }
    
    Info << endl;
    label patchID(0);
    word patchName("lowerWall");
    patchID =mesh.boundaryMesh().findPatchID(patchName);

    Info << "Found patch: " << patchName <<" at the ID " << patchID << nl <<endl;

    scalar surfaceArea = 0;
    forAll(mesh.boundary()[patchID].Sf(),faceI)
        surfaceArea += mag(mesh.boundary()[patchID].Sf()[faceI]);
    Info << "Surface Area of " << patchName << " = " <<surfaceArea <<endl;
*/
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
