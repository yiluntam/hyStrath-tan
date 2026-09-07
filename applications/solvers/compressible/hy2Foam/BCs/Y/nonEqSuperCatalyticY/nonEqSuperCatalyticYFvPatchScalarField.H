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
    Foam::nonEqSuperCatalyticYFvPatchScalarField

Description
    Supercatalytic wall boundary condition for a species mass fraction:
    the wall is assumed to be in instantaneous chemical equilibrium, so
    that the gas composition at the wall is prescribed (Dirichlet).
    In practice a fully recombined composition is prescribed (e.g.
    Y_N2 = 1, Y_N = 0 in a N2/N flow), which is the gamma -> infinity,
    supply-limited limit: the wall diffusive fluxes supply whatever the
    equilibrium requires.

    超催化壁面组分边界条件：假设壁面化学瞬时平衡，直接给定壁面组分
    （Dirichlet）。实际使用中给定完全复合组分（如 N2/N 流动中
    Y_N2 = 1、Y_N = 0），对应 γ → ∞ 的供给受限极限。

    This is distinct from the finite-rate fully catalytic limit
    (nonEqScottCatalyticY with gamma = 1), which is the kinetic ceiling
    gamma*c_bar/4 of the collision rate. The two merge in the
    continuum-thick-boundary-layer regime and separate with increasing
    rarefaction (an experimental validation discriminator).

    注意与有限速率完全催化（nonEqScottCatalyticY 且 γ = 1，碰撞速率的
    动力学上限 γ·c̄/4）相区别：二者在连续介质厚边界层区重合，随稀薄度
    增大而分离（可作为实验判别依据）。

    Because all species are prescribed, the impermeability of the wall is
    enforced automatically by the Stefan correction of the solver's
    diffusion model (sum of the Stefan-corrected wall fluxes is zero when
    the prescribed mass fractions sum to one).

    所有组分均给定，壁面不穿透性由求解器扩散模型的 Stefan 修正自动
    保证（给定质量分数之和为 1 时，Stefan 修正后的壁面通量之和为零）。

    'compositionMode' is a reserved keyword for a future equilibrium
    composition mode; only 'fixed' (the default, prescribed value) is
    currently available.

Usage
    \table
        Property        | Description            | Required | Default
        value           | 给定壁面质量分数        | yes      | -
        compositionMode | 预留关键字（仅 fixed） | no       | fixed
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nonEqSuperCatalyticY;
        value           uniform 1;   // e.g. N2 on a N2/N flow
    }
    \endverbatim

See also
    Foam::nonEqScottCatalyticYFvPatchScalarField
    Foam::nonEqScottCatalyticProductYFvPatchScalarField

SourceFiles
    nonEqSuperCatalyticYFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nonEqSuperCatalyticYFvPatchScalarField_H
#define nonEqSuperCatalyticYFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class nonEqSuperCatalyticYFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nonEqSuperCatalyticYFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("nonEqSuperCatalyticY");


    // Constructors

        //- Construct from patch and internal field
        nonEqSuperCatalyticYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nonEqSuperCatalyticYFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given nonEqSuperCatalyticYFvPatchScalarField
        //  onto a new patch
        nonEqSuperCatalyticYFvPatchScalarField
        (
            const nonEqSuperCatalyticYFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nonEqSuperCatalyticYFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nonEqSuperCatalyticYFvPatchScalarField
        (
            const nonEqSuperCatalyticYFvPatchScalarField&,
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
                new nonEqSuperCatalyticYFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
