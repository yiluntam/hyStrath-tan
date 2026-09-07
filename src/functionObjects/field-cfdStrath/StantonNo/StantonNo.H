/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
    Foam::functionObjects::StantonNo

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs the wall Stanton number for wall patches.

    The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

    Definition (classical hypersonic form, Fay-Riddell):

        St = q_w / (rho_inf*U_inf*(h_r - h_w))

    where the recovery enthalpy difference is approximated with the
    trans-rotational heat capacity:

        h_r - h_w = cp_inf*(T_r - T_w)
        T_r       = T_inf + r*U_inf^2/(2*cp_inf)   (recovery temperature)

    Reference quantities (rho_inf, U_inf, T_inf, cp_inf) are taken from
    the first face of the inflow patch; the recovery factor r defaults
    to sqrt(Pr_inf) (laminar boundary layer) with
    Pr_inf = mu_inf*cp_inf/kappa_tr,inf, or can be prescribed directly.

    The wall temperature T_w is read face-by-face from thermo.T(), so it
    follows whatever wall thermal condition is applied (isothermal,
    Smoluchowski jump, ...).

    Note: the field must be re-read from the start-time directory (see
    execute()); restart the case (startFrom latestTime) so that the
    wallHeatFlux file contains physical values before meaningful St can
    be produced.

Usage
    Example of function object specification to calculate the StantonNo:
    \verbatim
    StantonNo1
    {
        type        StantonNo;
        libs        ("libstrathFieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property        | Description                          | Required | Default value
        type            | Type name: StantonNo                 | yes      |
        wallHeatFlux    | Name of wallHeatFlux field           | no       | wallHeatFlux
        inflowPatchName | Name of inflow patch                 | no       | inlet
        recoveryFactor  | Recovery factor r; < 0 = auto (sqrt(Pr_inf), laminar) | no | -1
        result          | Name of Stanton field                | no       | \<function name\>
        log             | Log to standard output               | no       | yes
    \endtable


See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    StantonNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_StantonNo_H
#define functionObjects_StantonNo_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class StantonNo Declaration
\*---------------------------------------------------------------------------*/

class StantonNo
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        //- Name of wallHeatFlux field, default is "wallHeatFlux"
        word wallHeatFlux_;

        //- Name of inflow patch
        word inflowPatchName_;

        //- Recovery factor r; < 0 = auto (sqrt(Pr_inf), laminar BL)
        scalar recoveryFactor_;

        //- header of wallHeatFlux field
        IOobject wallHeatFluxHeader_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        StantonNo(const StantonNo&);

        //- Disallow default bitwise assignment
        void operator=(const StantonNo&);


public:

    //- Runtime type information
    TypeName("StantonNo");


    // Constructors

        //- Construct from Time and dictionary
        StantonNo
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~StantonNo();


    // Member Functions

        //- Read the StantonNo data
        virtual bool read(const dictionary&);

        //- Calculate the StantonNo field
        virtual bool execute();

        //- Write the StantonNo field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
