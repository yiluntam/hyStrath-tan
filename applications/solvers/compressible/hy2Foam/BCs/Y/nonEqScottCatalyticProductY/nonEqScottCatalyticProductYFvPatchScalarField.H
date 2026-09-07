/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\  /    A nd           | Copyright (C) 2016-2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    Foam::nonEqScottCatalyticProductYFvPatchScalarField

Description
    Catalytic recombination-product wall boundary condition, to be used on
    the molecular species produced by the surface recombination (e.g. N2
    for N recombination) on patches where the atomic species uses
    nonEqScottCatalyticY.

    复合产物组分的催化壁面边界条件，与原子组分的 nonEqScottCatalyticY
    配套用于同一壁面（例如 N 复合时用于 N2）。

    The recombined mass re-enters the gas phase as a diffusive influx of
    the product species, omega taken from the atom BC:
        snGrad(Y_product) = omega/rhoD
    i.e. a mixed condition with valueFraction = 0, refValue = 0,
    refGrad = omega/rhoD. Together with the atom Robin condition this
    closes the wall mass balance (sum of the wall diffusive fluxes is
    zero).

    复合质量以产物组分的扩散入流通量回到气相，ω 取自原子 BC：
    snGrad(Y_产物) = ω/ρD（valueFraction = 0 的 mixed 形式）。
    与原子 Robin 条件联立闭合壁面质量平衡（壁面扩散通量之和为零）。

Usage
    \table
        Property   | Description                       | Required | Default
        atomField  | 原子组分场名（须为 nonEqScottCatalyticY） | yes | -
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nonEqScottCatalyticProductY;
        atomField       N;
        value           uniform 0.8;
    }
    \endverbatim

See also
    Foam::nonEqScottCatalyticYFvPatchScalarField

SourceFiles
    nonEqScottCatalyticProductYFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nonEqScottCatalyticProductYFvPatchScalarField_H
#define nonEqScottCatalyticProductYFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class nonEqScottCatalyticProductYFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nonEqScottCatalyticProductYFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of the recombining atomic species field 原子组分场名
        word atomFieldName_;


public:

    //- Runtime type information
    TypeName("nonEqScottCatalyticProductY");


    // Constructors

        //- Construct from patch and internal field
        nonEqScottCatalyticProductYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonEqScottCatalyticProductYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonEqScottCatalyticProductYFvPatchScalarField
        //  onto a new patch
        nonEqScottCatalyticProductYFvPatchScalarField
        (
            const nonEqScottCatalyticProductYFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nonEqScottCatalyticProductYFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonEqScottCatalyticProductYFvPatchScalarField
        (
            const nonEqScottCatalyticProductYFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nonEqScottCatalyticProductYFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap(const fvPatchField<scalar>&, const labelList&);


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
