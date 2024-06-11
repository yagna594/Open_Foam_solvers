/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::doubleVelocityFvPatchVectorField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal vector boundary condition
    by its magnitude.

Usage
    \table
        Property     | Description                          | Required | Default
        refValue     | reference value                      | yes |
        ramp         | time-based ramping                   | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceNormalFixedValue;
        refValue        uniform -10;           // 10 INTO the domain
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors
    ControlDict Specification:
    -  libs ("libdoubleVelocityBC.so");

See also
    Foam::fixedValueFvPatchField

SourceFiles
    doubleVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef doubleVelocityFvPatchVectorField_H
#define doubleVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class doubleVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class doubleVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        scalarField refValue_;

        //- Optional time ramping
        autoPtr<Function1<scalar>> ramp_;


public:

    //- Runtime type information
    TypeName("libdoubleVelocityBC");


    // Constructors

        //- Construct from patch and internal field
        doubleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        doubleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  doubleVelocityFvPatchVectorField
        //  onto a new patch
        doubleVelocityFvPatchVectorField
        (
            const doubleVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        doubleVelocityFvPatchVectorField
        (
            const doubleVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new doubleVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        doubleVelocityFvPatchVectorField
        (
            const doubleVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new doubleVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::ftlLogProfileFvPatchVectorField

Group
    grpGenericBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a surface-normal vector boundary condition
    for an atmoaspheric boundary layer profile inlet. This is calculated from a 
    refernce height and reference velocity.

Usage
    \table
        Property     | Description                          | Required | Default
        refHeight    | reference Height                     | yes |
        refVelocity  | reference Velocity                   | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceNormalFixedValue;
        refHeight       100;           
        refVelocity     10;       
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors

    ControlDict Specification:
    -  libs(*libftlLogProfile.so*);

See also
    Foam::fixedValueFvPatchField

SourceFiles
    ftlLogProfileFvPatchVectorField.C

\*---------------------------------------------------------------------------

#ifndef doubleVelocityFvPatchVectorField_H
#define doubleVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

namespace Foam
{
    class doubleVelocityFvPatchVectorField : public fixedValueFvPatchVectorField
    {
    public:
        // Constructors
        doubleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        doubleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        // Update coefficients
        virtual void updateCoeffs();

        // Clone functions
        virtual tmp<fvPatchVectorField> clone() const;
        virtual tmp<fvPatchVectorField> clone(const DimensionedField<vector, volMesh>& iF) const;
    };
}

#endif

*/