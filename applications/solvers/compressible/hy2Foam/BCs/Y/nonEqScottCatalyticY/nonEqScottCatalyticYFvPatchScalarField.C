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

#include "nonEqScottCatalyticYFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"
#include "constants.H"

#include "nonEqSmoluchowskiJumpTFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonEqScottCatalyticYFvPatchScalarField::
nonEqScottCatalyticYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    rhoName_("rho"),
    TtName_("Tt"),
    catalyticEfficiency_(0.0),
    molarMass_(1.0),
    Twall_(p.size(), 0.0),
    TwallFromDict_(false),
    warned_(false),
    kc_(p.size(), 0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::nonEqScottCatalyticYFvPatchScalarField::
nonEqScottCatalyticYFvPatchScalarField
(
    const nonEqScottCatalyticYFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    rhoName_(ptf.rhoName_),
    TtName_(ptf.TtName_),
    catalyticEfficiency_(ptf.catalyticEfficiency_),
    molarMass_(ptf.molarMass_),
    Twall_(ptf.Twall_, mapper),
    TwallFromDict_(ptf.TwallFromDict_),
    warned_(false),
    kc_(ptf.kc_, mapper)
{}


Foam::nonEqScottCatalyticYFvPatchScalarField::
nonEqScottCatalyticYFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    TtName_(dict.lookupOrDefault<word>("Tt", "Tt")),
    catalyticEfficiency_(readScalar(dict.lookup("catalyticEfficiency"))),
    molarMass_(readScalar(dict.lookup("molarMass"))),
    Twall_(p.size(), 0.0),
    TwallFromDict_(false),
    warned_(false),
    kc_(p.size(), 0.0)
{
    if
    (
        catalyticEfficiency_ < 0.0
     || catalyticEfficiency_ > 1.0
    )
    {
        FatalIOErrorIn
        (
            "nonEqScottCatalyticYFvPatchScalarField::"
            "nonEqScottCatalyticYFvPatchScalarField"
            "("
            "    const fvPatch&,"
            "    const DimensionedField<scalar, volMesh>&,"
            "    const dictionary&"
            ")",
            dict
        )   << "unphysical catalyticEfficiency specified"
            << " (0 <= gamma <= 1)" << endl
            << exit(FatalIOError);
    }

    if (molarMass_ <= 0.0)
    {
        FatalIOErrorIn
        (
            "nonEqScottCatalyticYFvPatchScalarField::"
            "nonEqScottCatalyticYFvPatchScalarField"
            "("
            "    const fvPatch&,"
            "    const DimensionedField<scalar, volMesh>&,"
            "    const dictionary&"
            ")",
            dict
        )   << "molarMass must be positive [kg/kmol]" << endl
            << exit(FatalIOError);
    }

    // Optional prescribed wall temperature; if absent, the wall temperature
    // is resolved in updateCoeffs() from the Tt boundary condition so that
    // the catalytic wall and the temperature-jump wall share the same
    // surface temperature (single source).
    // 可选的给定壁温；缺省时在 updateCoeffs() 中从 Tt 边界条件解析，
    // 使催化壁面与温度跳跃模型共用同一壁面温度（同源）。
    if (dict.found("Twall"))
    {
        Twall_ = scalarField("Twall", dict, p.size());
        TwallFromDict_ = true;
    }

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


Foam::nonEqScottCatalyticYFvPatchScalarField::
nonEqScottCatalyticYFvPatchScalarField
(
    const nonEqScottCatalyticYFvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    rhoName_(ptpsf.rhoName_),
    TtName_(ptpsf.TtName_),
    catalyticEfficiency_(ptpsf.catalyticEfficiency_),
    molarMass_(ptpsf.molarMass_),
    Twall_(ptpsf.Twall_),
    TwallFromDict_(ptpsf.TwallFromDict_),
    warned_(false),
    kc_(ptpsf.kc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::nonEqScottCatalyticYFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    Twall_.autoMap(m);
    kc_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::nonEqScottCatalyticYFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const nonEqScottCatalyticYFvPatchScalarField& srptf =
        refCast<const nonEqScottCatalyticYFvPatchScalarField>(ptf);

    Twall_.rmap(srptf.Twall_, addr);
    kc_.rmap(srptf.kc_, addr);
}


// Update the coefficients associated with the patch field
void Foam::nonEqScottCatalyticYFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Resolve the wall temperature: dictionary entry, then the Tt BC's
    // wall temperature (temperature-jump condition), then a fixedValue Tt,
    // else the gas-side Tt patch value with a warning.
    // 壁温解析：字典项 > 温度跳跃条件的壁温 > fixedValue 的 Tt >
    // （带警告的）气侧 Tt 面值。
    scalarField Tw(patch().size(), 0.0);

    if (TwallFromDict_)
    {
        Tw = Twall_;
    }
    else
    {
        const volScalarField& Tt = db().lookupObject<volScalarField>(TtName_);

        const fvPatchScalarField& pTt =
            Tt.boundaryField()[patch().index()];

        if (isA<nonEqSmoluchowskiJumpTFvPatchScalarField>(pTt))
        {
            Tw =
                refCast<const nonEqSmoluchowskiJumpTFvPatchScalarField>(pTt)
                .Twall();
        }
        else if (isA<fixedValueFvPatchScalarField>(pTt))
        {
            Tw = pTt;
        }
        else
        {
            Tw = pTt;

            if (!warned_)
            {
                WarningIn
                (
                    "nonEqScottCatalyticYFvPatchScalarField::updateCoeffs()"
                )   << "Catalytic wall patch '" << patch().name()
                    << "': the Tt boundary condition is neither"
                    << " nonEqSmoluchowskiJumpT nor fixedValue."
                    << " The gas-side Tt patch value is used as the wall"
                    << " surface temperature. Prescribe 'Twall' in the"
                    << " catalytic BC dictionary or use a temperature-jump /"
                    << " fixedValue Tt condition for a consistent surface"
                    << " temperature." << endl;
                warned_ = true;
            }
        }
    }

    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    // Binary diffusion coefficient of this atomic species: the same field
    // the Y-equation laplacian uses ("rhoD_<species>").
    // 与 Y 方程 laplacian 相同的二组分扩散系数场（rhoD_<species>）。
    const fvPatchScalarField& prhoD =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "rhoD_" + internalField().name()
        );

    // Specific gas constant of the atomic species [J/(kg K)]:
    // RR = universal gas constant [J/(kmol K)], W = molar mass [kg/kmol]
    const scalar RR = 1000.0*constant::physicoChemical::R.value();

    // Surface reaction-rate coefficient, k_c = gamma*sqrt(Rs*Tw/(2*pi))
    // 表面反应速率系数 k_c = γ·√(R_s·T_w/(2π)) [m/s]
    kc_ =
        catalyticEfficiency_
       *sqrt(RR/(2.0*constant::mathematical::pi*molarMass_)*Tw);

    // Wall mass balance of the atoms (Robin condition):
    //     -rhoD*snGrad(Y) = k_c*rho*Y_w
    // expressed in mixed form:
    //     valueFraction = k_c*rho / (k_c*rho + rhoD*deltaCoeffs)
    // gamma = 0 gives valueFraction = 0, i.e. exactly zeroGradient.
    // 原子壁面质量平衡（Robin 条件）的 mixed 表达；γ = 0 时退化为
    // 精确的 zeroGradient（非催化逐位回归不变）。
    valueFraction() =
        kc_*prho/max(kc_*prho + prhoD*patch().deltaCoeffs(), SMALL);
    refValue() = 0.0;
    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::nonEqScottCatalyticYFvPatchScalarField::massFlux() const
{
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    // Wall mass fraction from the mixed formulation with the current
    // coefficients (identical to what evaluate() produces)
    // 用当前系数按 mixed 公式重构壁面质量分数（与 evaluate() 一致）
    scalarField Yw
    (
        valueFraction()*refValue()
      + (1.0 - valueFraction())
       *(patchInternalField() + refGrad()/patch().deltaCoeffs())
    );

    // omega = k_c*rho*Y_w [kg/(m^2 s)], positive: atoms recombine at the
    // wall, i.e. mass leaves the gas phase
    // 复合质量通量，正值表示原子在壁面复合（质量离开气相）
    return kc_*prho*Yw;
}


// Write
void Foam::nonEqScottCatalyticYFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "Tt", "Tt", TtName_);

    os.writeKeyword("catalyticEfficiency")
        << catalyticEfficiency_ << token::END_STATEMENT << nl;
    os.writeKeyword("molarMass")
        << molarMass_ << token::END_STATEMENT << nl;

    if (TwallFromDict_)
    {
        Twall_.writeEntry("Twall", os);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nonEqScottCatalyticYFvPatchScalarField
    );
}

// ************************************************************************* //
