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

#include "doubleVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::doubleVelocityFvPatchVectorField::
doubleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(p.size()),
    ramp_(nullptr)
{}


Foam::doubleVelocityFvPatchVectorField::
doubleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    refValue_("refValue", dict, p.size()),
    ramp_(Function1<scalar>::NewIfPresent("ramp", dict))
{
    tmp<vectorField> tvalues(refValue_*patch().nf());

    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(this->db().time().timeOutputValue());
    }

    fvPatchVectorField::operator=(tvalues);
}


Foam::doubleVelocityFvPatchVectorField::
doubleVelocityFvPatchVectorField
(
    const doubleVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(ptf.refValue_, mapper, pTraits<scalar>::zero),
    ramp_(ptf.ramp_.clone())
{
    // Note: refValue_ will have default value of 0 for unmapped faces. This
    // can temporarily happen during e.g. redistributePar. We should not
    // access ptf.patch() instead since redistributePar has destroyed this
    // at the time of mapping.

    tmp<vectorField> tvalues(refValue_*patch().nf());

    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(this->db().time().timeOutputValue());
    }

    fvPatchVectorField::operator=(tvalues);
}


Foam::doubleVelocityFvPatchVectorField::
doubleVelocityFvPatchVectorField
(
    const doubleVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    refValue_(ptf.refValue_),
    ramp_(ptf.ramp_.clone())
{}


Foam::doubleVelocityFvPatchVectorField::
doubleVelocityFvPatchVectorField
(
    const doubleVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    refValue_(ptf.refValue_),
    ramp_(ptf.ramp_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::doubleVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMap(mapper);
    refValue_.autoMap(mapper);
}


void Foam::doubleVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const doubleVelocityFvPatchVectorField& tiptf =
        refCast<const doubleVelocityFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::doubleVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<vectorField> tvalues(refValue_*patch().nf());
    if (ramp_)
    {
        tvalues.ref() *= ramp_->value(this->db().time().timeOutputValue());
    }
    
    fvPatchVectorField::operator=(2*tvalues);
    fvPatchVectorField::updateCoeffs();
}


void Foam::doubleVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    refValue_.writeEntry("refValue", os);
    if (ramp_)
    {
        ramp_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        doubleVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
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

\*---------------------------------------------------------------------------

#include "doubleVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// Constructors
Foam::doubleVelocityFvPatchVectorField::doubleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}

Foam::doubleVelocityFvPatchVectorField::doubleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}

vectorField tvalues = patch().nf();
    forAll(tvalues,faceI)
    {
        tvalues[faceI] =  2 * refValue_
        //-1*(u_star/K)*log((zHeight[faceI][2]+z0)/z0);
    }
// Update coefficients
void Foam::doubleVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Double the specified velocity value
    operator==(2.0 * refValue());

    fixedValueFvPatchVectorField::updateCoeffs();
}

// Clone functions
tmp<fvPatchVectorField> Foam::doubleVelocityFvPatchVectorField::clone() const
{
    return tmp<fvPatchVectorField>
    (
        new doubleVelocityFvPatchVectorField(*this)
    );
}

tmp<fvPatchVectorField> Foam::doubleVelocityFvPatchVectorField::clone
(
    const DimensionedField<vector, volMesh>& iF
) const
{
    return tmp<fvPatchVectorField>
    (
        new doubleVelocityFvPatchVectorField(patch(), iF)
    );
}

// Register the boundary condition
makePatchTypeField
(
    fvPatchVectorField,
    doubleVelocityFvPatchVectorField
);
*/