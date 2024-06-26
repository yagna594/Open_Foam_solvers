/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "ftlLogProfileFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ftlLogProfileFvPatchVectorField::
ftlLogProfileFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refHeight_(0.),
    refVelocity_(0.)
{}


Foam::ftlLogProfileFvPatchVectorField::
ftlLogProfileFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    refHeight_(0.),
    refVelocity_(0.)
{
    dict.lookup("refHeight") >> refHeight_;
    dict.lookup("refVelocity") >> refVelocity_;

    updateCoeffs();
}


Foam::ftlLogProfileFvPatchVectorField::
ftlLogProfileFvPatchVectorField
(
    const ftlLogProfileFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    refHeight_(ptf.refHeight_),
    refVelocity_(ptf.refVelocity_)
{
    updateCoeffs();
}


Foam::ftlLogProfileFvPatchVectorField::
ftlLogProfileFvPatchVectorField
(
    const ftlLogProfileFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    refHeight_(ptf.refHeight_),
    refVelocity_(ptf.refVelocity_)
{}


Foam::ftlLogProfileFvPatchVectorField::
ftlLogProfileFvPatchVectorField
(
    const ftlLogProfileFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    refHeight_(ptf.refHeight_),
    refVelocity_(ptf.refVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ftlLogProfileFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMap(mapper);
}


void Foam::ftlLogProfileFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const ftlLogProfileFvPatchVectorField& tiptf =
        refCast<const ftlLogProfileFvPatchVectorField>(ptf);
}


void Foam::ftlLogProfileFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Get the Z height of the centroid coordinates of the cells within the patch
    const vectorField& zHeight = patch().Cf();

    vectorField tvalues = patch().nf();
    scalar K = 0.40;
    scalar z0 = 0.3;

    scalar u_star = K * (refVelocity_/log((refHeight_+z0)/z0));
    forAll(tvalues,faceI)
    {
        tvalues[faceI] *= -1*(u_star/K)*log((zHeight[faceI][2]+z0)/z0);
    }
    //Assign the computed values to the boundary condition
    fvPatchVectorField::operator=(tvalues);
    //Call the base class method to update coefficients
    fvPatchVectorField::updateCoeffs();
}


void Foam::ftlLogProfileFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("refHeight") << refHeight_ <<token::END_STATEMENT<<nl;
    os.writeKeyword("refVelocity") << refVelocity_ <<token::END_STATEMENT<<nl;
    /*
    refValue_.writeEntry("refValue", os);
    if (ramp_)
    {
        ramp_->writeData(os);
    }*/
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        ftlLogProfileFvPatchVectorField
    );
}

// ************************************************************************* //
