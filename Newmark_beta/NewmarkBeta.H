/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

Class
    Foam::sixDoFSolvers::NewmarkBeta

Group
    grpSixDoFRigidBodySolvers

Description
    NewmarkBeta 2nd-order time-integrator for 6DoF solid-body motion.

    Reference:
    \verbatim
        NewmarkBeta, N. M. (1959).
        A method of computation for structural dynamics.
        Journal of the Engineering Mechanics Division, 85(3), 67-94.
    \endverbatim

    Example specification in dynamicMeshDict:
    \verbatim
    solver
    {
        type    NewmarkBeta;
        gamma   0.5;    // Velocity integration coefficient
        beta    0.25;   // Position integration coefficient
    }
    \endverbatim

SourceFiles
    NewmarkBeta.C

\*---------------------------------------------------------------------------*/

#ifndef NewmarkBeta_H
#define NewmarkBeta_H

#include "sixDoFSolver.H"
//#include "sixDoFRigidBodyMotionRestraint.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{

/*---------------------------------------------------------------------------*\
                           Class NewmarkBeta Declaration
\*---------------------------------------------------------------------------*/

class NewmarkBeta
:
    public sixDoFSolver
{
    // Private data

        //- Coefficient for velocity integration (default: 0.5)
        const scalar gamma_;

        //- Coefficient for position and orientation integration (default: 0.25)
        const scalar beta_;

       const scalar mass_;
       scalar damping_;
       scalar stiffness_;

public:

    //- Runtime type information
    TypeName("NewmarkBeta");


    // Constructors
       // NewmarkBeta
        //(
        //    const sixDoFRigidBodyMotion& sDoFRBM
        //);
 
        //- Construct from a dictionary and the body
        NewmarkBeta
        (
            const dictionary& dict,
            sixDoFRigidBodyMotion& body
        );

        //- Construct and return a clone
        virtual autoPtr<sixDoFSolver> clone() const
        {
            return autoPtr<sixDoFSolver>::NewFrom<NewmarkBeta>(dict_, body_);
        }


    //- Destructor
    virtual ~NewmarkBeta();


    // Member Functions

        //- Drag coefficient
        virtual void solve
        (
            bool firstIter,
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT,
            scalar deltaT0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sixDoFSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
