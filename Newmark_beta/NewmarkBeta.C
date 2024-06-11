/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/
//#include "fvCFD.H"
#include "NewmarkBeta.H"
#include "addToRunTimeSelectionTable.H"
//#include "createFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFSolvers
{
    defineTypeNameAndDebug(NewmarkBeta, 0);
    addToRunTimeSelectionTable(sixDoFSolver, NewmarkBeta, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Constructor implementation
/*Foam::NewmarkBeta::NewmarkBeta(const dictionary& dict):
{
    const dictionary& sixDoFRigidBodyMotionCoeffs = dict.subDict("sixDoFRigidBodyMotionCoeffs");

    // Extract mass, damping, and stiffness
    mass = sixDoFRigidBodyMotionCoeffs.lookupOrDefault<scalar>("mass", 1.0);
    
    if (sixDoFRigidBodyMotionCoeffs.found("restraints"))
    {
        const dictionary& restraints = sixDoFRigidBodyMotionCoeffs.subDict("restraints");
        const dictionary& verticalSpring = restraints.subDict("verticalSpring");

        damping = verticalSpring.lookupOrDefault<scalar>("damping", 0.0);
        stiffness = verticalSpring.lookupOrDefault<scalar>("stiffness", 0.0);
    }
    else
    {
        damping = 0.0;
        stiffness = 0.0;
    }
}*/
/*Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const sixDoFRigidBodyMotion& sDoFRBM
)
:
    mass_(sDoFRBM.mass_)
{}
*/
Foam::sixDoFSolvers::NewmarkBeta::NewmarkBeta
(
    const dictionary& dict,
    sixDoFRigidBodyMotion& body
)
:
    sixDoFSolver(dict, body),
    gamma_(dict.getOrDefault<scalar>("gamma", 0.5)),
    beta_
    (
        max
        (
            0.25*sqr(gamma_ + 0.5),
            dict.getOrDefault<scalar>("beta", 0.25)
        )
    ),
    mass_(body.mass()),
    damping_(0.0),
    stiffness_(0.0)
{
    if (dict.found("sixDoFRigidBodyMotionCoeffs"))
    {
        const dictionary& sixDoFRigidBodyMotionCoeffs = dict.subDict("sixDoFRigidBodyMotionCoeffs");

        if (sixDoFRigidBodyMotionCoeffs.found("restraints"))
        {
            const dictionary& restraints = sixDoFRigidBodyMotionCoeffs.subDict("restraints");

            if (restraints.found("verticalSpring"))
            {
                const dictionary& verticalSpring = restraints.subDict("verticalSpring");

                damping_ = verticalSpring.lookupOrDefault<scalar>("damping", 0.0);
                stiffness_ = verticalSpring.lookupOrDefault<scalar>("stiffness", 0.0);
            }
        }
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFSolvers::NewmarkBeta::~NewmarkBeta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*//Mass
scalar mass;
dynamicMeshDict.lookup("mass") >> mass;
//Spring stiffness
scalar stiffness;
dynamicMeshDict.lookup("stiffness") >> stiffness;
//Reference pressure value
scalar damping;
dynamicMeshDict.lookup("damping") >> damping;*/

void Foam::sixDoFSolvers::NewmarkBeta::solve
(
    bool firstIter,
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT,
    scalar deltaT0
)
{
    
    /*// Update the linear acceleration and torque
    updateAcceleration(fGlobal, tauGlobal);

    // Define NewmarkBeta-beta parameters
    scalar gamma = 0.5;  // NewmarkBeta-gamma
    scalar beta = 0.25;  // NewmarkBeta-beta

   // a = acceleration(n+1) a0 = acceleration(n)
    // Correct linear velocity using NewmarkBeta-beta method
    // velocity(n+1) = velocity(n) + (1-gamma) deltat acceleration + gamma.deltat.acceleration(n+1)
    v() = v0() + deltaT * (1.0 - gamma) * a0() + deltaT * gamma * a();

    // Correct angular momentum using NewmarkBeta-beta method
    pi() = pi0() + deltaT * (1.0 - gamma) * tau0() + deltaT * gamma * tau();
    
    // Correct position using NewmarkBeta-beta method    
    // displacement(n+1) = displacement(n) +deltat velocity(n) +deltat^2.0.5.((1-2beta)acceleration(n)+2.betta.acceleration(n+1))   

    centreOfRotation() = centreOfRotation0() + deltaT * v0() + 0.5 * sqr(deltaT) * ((1.0 - 2.0 * beta) * a0() + 2.0 * beta * a());

    // Correct orientation
    vector piDeltaT = deltaT * pi0() + 0.5 * sqr(deltaT) * ((1.0 - 2.0 * beta) * tau0() + 2.0 * beta * tau());
    Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
    Q() = Qpi.first();
    */
    //scalar mass = mass_(dict.get<scalar>("mass"));
    //scalar damping = damping_(dict.get<scalar>("damping"));
    //scalar stiffness = stiffness_(dict.get<scalar>("stiffness"));
    scalar gamma = 0.5;  // NewmarkBeta-gamma
    scalar beta = 0.25;  // NewmarkBeta-beta

    // Predictor step
    vector v_pred = v() + deltaT * (1.0 - gamma) * a0();
    vector pi_pred = pi() + deltaT * (1.0 - gamma) * tau0();
    vector centreOfRotation_pred = centreOfRotation() + deltaT * v0() +
                                   0.5 * sqr(deltaT) * ((1.0 - 2.0 * beta) * a0());

    // Update acceleration and torque using the predicted values
    scalar denom = mass_ + damping_ * gamma_ * deltaT + stiffness_ * beta_ * deltaT * deltaT;
    if (denom == 0)
    {
        FatalErrorInFunction
            << "Denominator is zero, which will cause division by zero." << nl
            << "mass: " << mass_ << nl
            << "damping: " << damping_ << nl
            << "stiffness: " << stiffness_ << nl
            << "gamma: " << gamma_ << nl
            << "beta: " << beta_ << nl
            << "deltaT: " << deltaT << exit(Foam::FatalError);
    }
    a() = (fGlobal - damping_ * v_pred - stiffness_ * centreOfRotation_pred) / denom;
    tau() = tauGlobal; 
    // Assuming tauGlobal is updated elsewhere in the code
    // Update acceleration and torque using the predicted values
    // ddu = (f - c*du_pred - k*u_pred) / (m + c * gamma * deltaT + k * beta * deltaT^2)
    //a() = (fGlobal - damping_ * v_pred - stiffness_ * centreOfRotation_pred) / (mass_ + damping_ * gamma * deltaT + stiffness_ * beta * sqr(deltaT));
    
    a() = tConstraints() & a();
    // Correct linear velocity using NewmarkBeta-beta method
    // velocity(n+1) = velocity(n) + (1-gamma) deltat acceleration + gamma.deltat.acceleration(n+1)
    v() = tConstraints() & (v_pred + deltaT * gamma * a());

    // Correct angular momentum using NewmarkBeta-beta method
    pi() = rConstraints() & (pi_pred + deltaT * gamma * tau());

    // Correct position using NewmarkBeta-beta method    
    // displacement(n+1) = displacement(n) + deltat velocity(n) + deltat^2 * 0.5 * ((1 - 2 * beta) * acceleration(n) + 2 * beta * acceleration(n+1))
    centreOfRotation() = centreOfRotation_pred + (tConstraints() &  (0.5 * sqr(deltaT) * (2.0 * beta * a()))); 
    // Correct orientation
    vector piDeltaT = rConstraints() &  (deltaT * pi0() + 0.5 * sqr(deltaT) * ((1.0 - 2.0 * beta) * tau0() + 2.0 * beta * tau()));
    Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
    Q() = Qpi.first();
    /****************************************************************************************************************************/
/*    // Correct linear velocity
    v() =
        tConstraints()
      & (v0() + aDamp()*deltaT*(gamma_*a() + (1 - gamma_)*a0()));

    // Correct angular momentum
    pi() =
        rConstraints()
      & (pi0() + aDamp()*deltaT*(gamma_*tau() + (1 - gamma_)*tau0()));

    // Correct position
    centreOfRotation() =
        centreOfRotation0()
      + (
            tConstraints()
          & (
                deltaT*v0()
              + aDamp()*sqr(deltaT)*(beta_*a() + (0.5 - beta_)*a0())
            )
        );

    // Correct orientation
    vector piDeltaT =
        rConstraints()
      & (
            deltaT*pi0()
          + aDamp()*sqr(deltaT)*(beta_*tau() + (0.5 - beta_)*tau0())
        );
    Tuple2<tensor, vector> Qpi = rotate(Q0(), piDeltaT, 1);
    Q() = Qpi.first();*/
}

// ************************************************************************* //
