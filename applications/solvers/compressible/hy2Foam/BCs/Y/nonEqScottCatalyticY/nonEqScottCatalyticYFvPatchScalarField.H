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
    Foam::nonEqScottCatalyticYFvPatchScalarField

Description
    Finite-rate catalytic wall boundary condition for an atomic species
    (Scott model, gamma-collision form).

    原子组分的有限速率催化壁面边界条件（Scott 模型，γ 碰撞形式）。

    The surface recombination consumes atoms at the rate
        omega = gamma*rho_w*Y_w*sqrt(Rs*Tw/(2*pi)) = k_c*rho_w*Y_w,
    where k_c [m/s] is the surface reaction-rate coefficient and Rs = RR/W
    the specific gas constant of the atomic species. Together with the
    no-penetration wall mass balance this gives the Robin condition
        -rhoD*snGrad(Y) = k_c*rho*Y_w,
    implemented as a mixed (Robin) boundary condition with
        valueFraction = k_c*rho / (k_c*rho + rhoD*deltaCoeffs),
        refValue = 0, refGrad = 0.
    The recombination product species must use
    nonEqScottCatalyticProductY on the same patch so that the recombined
    mass re-enters the gas phase (wall mass balance).

    表面复合以 omega = γ·ρ_w·Y_w·√(R_s·T_w/(2π)) 的速率消耗原子；
    与壁面不穿透质量平衡联立得 Robin 条件 −ρD·∂Y/∂n = k_c·ρ·Y_w，
    以 mixed 形式实现。同一壁面上复合产物组分须配套使用
    nonEqScottCatalyticProductY，使复合质量回到气相（壁面质量闭合）。

    The wall temperature is taken, by priority:
    壁面温度按优先级取用：
      1. the 'Twall' entry of this BC dictionary, if present
         （本字典中的 Twall 项，若存在）
      2. the wall temperature of the Tt boundary condition if it is a
         nonEqSmoluchowskiJumpT (temperature-jump) condition
         （Tt 为温度跳跃条件时的壁温，与跳跃模型同源）
      3. the Tt boundary value if Tt is a fixedValue condition
         （Tt 为 fixedValue 时的壁面值）
      4. otherwise the gas-side Tt patch value, with a warning
         （否则用气侧 Tt 面值并给出警告）

    gamma = 0 gives valueFraction = 0, i.e. exactly zeroGradient
    (bitwise regression with non-catalytic settings).
    γ = 0 时 valueFraction = 0，退化为精确的 zeroGradient
    （非催化设置下逐位回归不变）。

Usage
    \table
        Property            | Description                | Required | Default
        catalyticEfficiency | 催化效率 γ (0 ≤ γ ≤ 1)     | yes      | -
        molarMass           | 原子摩尔质量 [kg/kmol]      | yes      | -
        Twall               | 壁面温度 [K]（均匀或逐面） | no       | from Tt BC
        rho                 | 密度场名                   | no       | rho
        Tt                  | 平动温度场名               | no       | Tt
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                nonEqScottCatalyticY;
        catalyticEfficiency 0.1;
        molarMass           14.0067;   // N
        value               uniform 0.2;
    }
    \endverbatim

See also
    Foam::nonEqScottCatalyticProductYFvPatchScalarField
    Foam::nonEqSuperCatalyticYFvPatchScalarField
    Foam::nonEqSmoluchowskiJumpTFvPatchScalarField

SourceFiles
    nonEqScottCatalyticYFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nonEqScottCatalyticYFvPatchScalarField_H
#define nonEqScottCatalyticYFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class nonEqScottCatalyticYFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nonEqScottCatalyticYFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Density field name, default = "rho" 密度场名
        word rhoName_;

        //- Translational-temperature field name, default = "Tt" 平动温度场名
        word TtName_;

        //- Catalytic efficiency gamma, 0 < gamma <= 1 催化效率
        scalar catalyticEfficiency_;

        //- Molar mass of the recombining atomic species [kg/kmol] 摩尔质量
        scalar molarMass_;

        //- Wall surface temperature prescribed in this dictionary, if any
        //  本字典给定的壁面温度（若存在）
        scalarField Twall_;

        //- Whether the wall temperature was read from the dictionary
        //  壁温是否来自字典
        bool TwallFromDict_;

        //- Warn-once flag for the gas-side-Twall fallback 警告一次标志
        mutable bool warned_;

        //- Surface reaction-rate coefficient
        //  k_c = gamma*sqrt(Rs*Tw/(2*pi)) [m/s] 表面反应速率系数
        scalarField kc_;


public:

    //- Runtime type information
    TypeName("nonEqScottCatalyticY");


    // Constructors

        //- Construct from patch and internal field
        nonEqScottCatalyticYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonEqScottCatalyticYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonEqScottCatalyticYFvPatchScalarField
        //  onto a new patch
        nonEqScottCatalyticYFvPatchScalarField
        (
            const nonEqScottCatalyticYFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nonEqScottCatalyticYFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonEqScottCatalyticYFvPatchScalarField
        (
            const nonEqScottCatalyticYFvPatchScalarField&,
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
                new nonEqScottCatalyticYFvPatchScalarField(*this, iF)
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


        // Access functions

            //- Return the catalytic recombination mass flux
            //  omega = k_c*rho*Y_w [kg/(m^2 s)], positive: atoms recombine
            //  at the wall (mass leaves the gas phase)
            //  催化复合质量通量，正值表示原子在壁面复合（质量离开气相）
            tmp<scalarField> massFlux() const;

            //- Return the surface reaction-rate coefficient [m/s]
            //  表面反应速率系数
            const scalarField& kc() const
            {
                return kc_;
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
