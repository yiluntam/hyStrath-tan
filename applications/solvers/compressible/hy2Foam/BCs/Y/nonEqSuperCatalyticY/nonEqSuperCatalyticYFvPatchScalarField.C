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

\*---------------------------------------------------------------------------*/

#include "nonEqSuperCatalyticYFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonEqSuperCatalyticYFvPatchScalarField::
nonEqSuperCatalyticYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::nonEqSuperCatalyticYFvPatchScalarField::
nonEqSuperCatalyticYFvPatchScalarField
(
    const nonEqSuperCatalyticYFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::nonEqSuperCatalyticYFvPatchScalarField::
nonEqSuperCatalyticYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
    // 'compositionMode' is reserved for a future equilibrium-composition
    // mode; only the prescribed-value mode is currently available.
    // compositionMode 为未来平衡组分模式预留的关键字。
    if (dict.found("compositionMode"))
    {
        const word compositionMode(dict.lookup("compositionMode"));

        if (compositionMode != "fixed")
        {
            FatalIOErrorIn
            (
                "nonEqSuperCatalyticYFvPatchScalarField::"
                "nonEqSuperCatalyticYFvPatchScalarField"
                "("
                "    const fvPatch&,"
                "    const DimensionedField<scalar, volMesh>&,"
                "    const dictionary&"
                ")",
                dict
            )   << "Unknown compositionMode '" << compositionMode
                << "'. Only 'fixed' (prescribed wall composition) is"
                << " currently available; 'equilibrium' is reserved for"
                << " future use." << endl
                << exit(FatalIOError);
        }
    }

    // Physical-range check of the prescribed wall mass fraction
    // 给定壁面质量分数的物理范围检查
    if (min(*this) < 0.0 || max(*this) > 1.0)
    {
        WarningIn
        (
            "nonEqSuperCatalyticYFvPatchScalarField::"
            "nonEqSuperCatalyticYFvPatchScalarField"
            "("
            "    const fvPatch&,"
            "    const DimensionedField<scalar, volMesh>&,"
            "    const dictionary&"
            ")"
        )   << "Prescribed wall mass fraction on patch '"
            << patch().name() << "' is outside [0, 1]: min = "
            << min(*this) << ", max = " << max(*this) << endl;
    }
}


Foam::nonEqSuperCatalyticYFvPatchScalarField::
nonEqSuperCatalyticYFvPatchScalarField
(
    const nonEqSuperCatalyticYFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Write
void Foam::nonEqSuperCatalyticYFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nonEqSuperCatalyticYFvPatchScalarField
    );
}

// ************************************************************************* //
