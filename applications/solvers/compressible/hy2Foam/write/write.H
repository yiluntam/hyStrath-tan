
// --- Output at write time
rarefactionParameters().write();
transport().write();
composition.write();

//- Write diffusion fluxes
if (not hyLight)
{
/*    volVectorField sumDiffusiveFluxes =
          speciesDiffusion().sumDiffusiveFluxes();*/
/*    sumDiffusiveFluxes.write();*/
/*    forAll(Y, speciei)*/
/*    {*/
/*        volVectorField Jc*/
/*        (*/
/*            IOobject*/
/*            (*/
/*                "J_" + composition.species()[speciei],*/
/*                runTime.timeName(),*/
/*                mesh,*/
/*                IOobject::NO_READ,*/
/*                IOobject::NO_WRITE*/
/*            ),*/
/*            speciesDiffusion().Jcorr(speciei)*/
/*        );*/
/*        Jc.write();*/
/*    }*/

/*    forAll(Y, speciei)*/
/*    {*/
/*        const volScalarField Di*/
/*        (*/
/*            IOobject*/
/*            (*/
/*                "D_" + composition.species()[speciei],*/
/*                runTime.timeName(),*/
/*                mesh,*/
/*                IOobject::NO_READ,*/
/*                IOobject::NO_WRITE*/
/*            ),*/
/*            speciesDiffusion().rhoD(speciei)*/
/*        );*/
/*        */
/*        Di.write();*/
/*    }*/
}

if (not downgradeSingleT)
{
    if (not hyLight)
    {
        thermo.Tov().write();

        if (noSolvedVibEq > 1)
        {
            thermo.Tv().write();
        }
    }
}

if (downgradeSingleTv)
{
    thermo.Tv().write();
}


if (writeWallHeatFlux)
{
    surfaceScalarField heatFlux_convtr =
        -fvc::interpolate(transport().kappatr())*fvc::snGrad(thermo.T());

    surfaceScalarField heatFlux_convve =
        -fvc::interpolate(transport().kappave())*fvc::snGrad(thermo.Tv());
        
    surfaceScalarField heatFlux_diff =
        speciesDiffusion().getDiffusiveWallHeatFlux();

    surfaceScalarField heatFlux = heatFlux_convtr + heatFlux_convve
        + heatFlux_diff;

    const surfaceScalarField::Boundary& patchHeatFlux =
        heatFlux.boundaryField();

    Info<< "\nIntegrated wall heat flux [W]" << endl;
    forAll(patchHeatFlux, patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            wallHeatFlux.boundaryFieldRef()[patchi] = patchHeatFlux[patchi];
            wallHeatFlux_tr.boundaryFieldRef()[patchi] = heatFlux_convtr.boundaryField()[patchi];
            wallHeatFlux_ve.boundaryFieldRef()[patchi] = heatFlux_convve.boundaryField()[patchi];

            const scalar int_tot_whf =
                gSum
                (
                    mesh.magSf().boundaryField()[patchi]
                   *patchHeatFlux[patchi]
                );
                
            scalar int_convtr_whf = 0.0;
            scalar int_convve_whf = 0.0;
            scalar int_diff_whf = 0.0;    

            if ( mag(int_tot_whf) > SMALL )
            {
                int_convtr_whf = 
                    gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                       *heatFlux_convtr.boundaryField()[patchi]
                    )/int_tot_whf*100.0;
                    
                int_convve_whf = 
                    gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                       *heatFlux_convve.boundaryField()[patchi]
                    )/int_tot_whf*100.0;
                    
                int_diff_whf =
                    gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                       *heatFlux_diff.boundaryField()[patchi]
                    )/int_tot_whf*100.0;  
            }
            
            Info<< "Patch " << patchi
                << " named " << patchHeatFlux[patchi].patch().name()
                << ": " << int_tot_whf << nl
                << "   convective trans-rot  (%): " << int_convtr_whf << nl
                << "   convective vibro-elec (%): " << int_convve_whf << nl
                << "   diffusive             (%): " << int_diff_whf << nl
                << endl;    
        }
    }
    Info<< endl;

    // Catalytic contribution to the wall heat flux (absolute-enthalpy
    // basis: includes the full recombination heat omega*Delta_h_recomb;
    // positive into the wall). Reported separately from the diffusive
    // item, which contains the same flux on the sensible/absolute basis
    // used by getDiffusiveWallHeatFlux().
    // 壁面热流的催化分量（绝对焓基准，含全部复合反应热；指向壁面
    // 为正），单独报告。
    if (catalyticWallPresent)
    {
        Info<< "Integrated catalytic wall heat flux [W]" << endl;

        forAll(catalyticPatchIDs, i)
        {
            const label patchi = catalyticPatchIDs[i];

            const scalar int_cat_whf =
                gSum
                (
                    mesh.magSf().boundaryField()[patchi]
                   *wallHeatFlux_cat.boundaryField()[patchi]
                );

            // Share of the total wall heat flux of the same patch
            // (Goulard-type diagnostic)
            // 同壁面总热流中的占比（Goulard 型诊断量）
            scalar int_cat_ratio_whf = 0.0;

            const scalar int_tot_patch_whf =
                gSum
                (
                    mesh.magSf().boundaryField()[patchi]
                   *wallHeatFlux.boundaryField()[patchi]
                );

            if (mag(int_tot_patch_whf) > SMALL)
            {
                int_cat_ratio_whf = int_cat_whf/int_tot_patch_whf*100.0;
            }

            Info<< "Patch " << patchi
                << " named " << mesh.boundaryMesh()[patchi].name()
                << ": " << int_cat_whf << nl
                << "   catalytic (% of total): " << int_cat_ratio_whf << nl
                << endl;
        }

        wallHeatFlux_cat.write();
    }

    wallHeatFlux.write();
    wallHeatFlux_tr.write();
    wallHeatFlux_ve.write();
}


if (reaction->lookupOrDefault<bool>("applyChangesAtWriteTime", false))
{
    return true;
}

if (reaction->lookupOrDefault<bool>("applyChangesAtWriteTimeAndWait", false))
{
    label restartIn = reaction->lookupOrDefault<label>("restartIn", 120);

    Info << "Waiting " << restartIn << " seconds before restart\n\n\n";

    Foam::sleep(restartIn);

    return true;
}
