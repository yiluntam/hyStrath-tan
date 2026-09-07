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

#include "nonEqScottCatalyticProductYFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "nonEqScottCatalyticYFvPatchScalarField.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonEqScottCatalyticProductYFvPatchScalarField::
nonEqScottCatalyticProductYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    atomFieldName_("atomFieldName")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::nonEqScottCatalyticProductYFvPatchScalarField::
nonEqScottCatalyticProductYFvPatchScalarField
(
    const nonEqScottCatalyticProductYFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    atomFieldName_(ptf.atomFieldName_)
{}


Foam::nonEqScottCatalyticProductYFvPatchScalarField::
nonEqScottCatalyticProductYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    atomFieldName_(dict.lookup("atomField"))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::nonEqScottCatalyticProductYFvPatchScalarField::
nonEqScottCatalyticProductYFvPatchScalarField
(
    const nonEqScottCatalyticProductYFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    atomFieldName_(ptpsf.atomFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::nonEqScottCatalyticProductYFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::nonEqScottCatalyticProductYFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void Foam::nonEqScottCatalyticProductYFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const volScalarField& atomField =
        db().lookupObject<volScalarField>(atomFieldName_);

    const fvPatchScalarField& patom =
        atomField.boundaryField()[patch().index()];

    if (!isA<nonEqScottCatalyticYFvPatchScalarField>(patom))
    {
        FatalErrorIn
        (
            "nonEqScottCatalyticProductYFvPatchScalarField::updateCoeffs()"
        )   << "Species '" << internalField().name()
            << "' on patch '" << patch().name()
            << "' uses nonEqScottCatalyticProductY, but the referenced"
            << " atomic species field '" << atomFieldName_
            << "' does not use nonEqScottCatalyticY on this patch."
            << endl
            << exit(FatalError);
    }

    // Make sure the atom BC coefficients are current. The call is
    // idempotent: updateCoeffs() is guarded by the updated() flag within
    // one boundary-update cycle, and the quantities it reads (wall
    // temperature, density, rhoD) do not change inside this cycle, so
    // re-evaluation later in the cycle returns identical coefficients.
    // 确保原子 BC 系数为最新。该调用幂等：updateCoeffs() 受 updated()
    // 标志保护，且其读取的量（壁温、密度、ρD）在本更新周期内不变。
    fvPatchScalarField& patomNC = const_cast<fvPatchScalarField&>(patom);
    refCast<nonEqScottCatalyticYFvPatchScalarField>(patomNC).updateCoeffs();

    // Recombination rate from the atom BC
    // 从原子 BC 取复合速率 ω [kg/(m^2 s)]
    const scalarField omega =
        refCast<const nonEqScottCatalyticYFvPatchScalarField>(patom)
        .massFlux();

    // Binary diffusion coefficient of this product species
    const fvPatchScalarField& prhoD =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "rhoD_" + internalField().name()
        );

    // The recombined mass re-enters the gas phase as a diffusive influx:
    //     J_product . n = -omega  (inward, n outward)
    //     -rhoD*snGrad(Y) = -omega -> snGrad(Y) = omega/rhoD
    // i.e. a pure gradient (valueFraction = 0) mixed condition.
    // If rhoD vanishes (trace species), the influx degenerates to a
    // zero gradient; the atom-side omega is self-limiting in that regime
    // (the atom wall value tends to zero).
    // 复合质量以扩散入流回到气相；ρD 趋于零时退化为零梯度（该极限下
    // 原子侧 ω 自身趋于零，自洽）。
    refGrad() = omega/max(prhoD, SMALL);
    refValue() = 0.0;
    valueFraction() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


// Write
void Foam::nonEqScottCatalyticProductYFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeKeyword("atomField")
        << atomFieldName_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nonEqScottCatalyticProductYFvPatchScalarField
    );
}

// ************************************************************************* //
