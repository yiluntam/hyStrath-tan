
// --- Solve energy, VISCOUS
if (!inviscid)
{
    if (downgradeSingleT)
    {
        // Single-temperature model: solve the TOTAL internal energy e with the
        // total thermal conductivity (upstream form). The trans-rotational
        // operator splitting (et + kappatr) is not applicable here: no separate
        // ev-equation is solved, so the VT-leakage rationale does not exist,
        // and kappatr would drop the kappave share of conduction. Solving e
        // also avoids the calculatedFvPatchField FatalError on et(=het), whose
        // single-T construction in rho2ReactionThermo carries no boundary types
        // (latent upstream bug, dormant upstream because the et boundary
        // gradient is never taken there).
        fvScalarMatrix eEqnViscous
        (
            fvm::ddt(rho, e) - fvc::ddt(rho, e)
          - fvm::laplacian(thermo.kappa()/thermo.CvMix(), e)
          ==
            fvOptions(rho, e)
        );

        // No separate ev-equation is solved, so the full species diffusion
        // enthalpy flux sum_j hs_j*J_j enters the energy budget here.
        if (speciesDiffusionName != "noSpeciesDiffusion")
        {
            eEqnViscous +=
                fvc::div
                (
                    speciesDiffusion().multiSpeciesHeatSource(),
                    "div(diffusiveFluxes)"
                );
        }

        // Catalytic wall: boundary completion of the species-diffusive
        // enthalpy flux (sensible share; the interior source term has
        // zero boundary values)
        // 催化壁面：组分扩散焓通量的边界完成项（显焓份额）
        if (catalyticWallPresent)
        {
            eEqnViscous +=
                fvc::div(catalyticFluxE(), "div(diffusiveFluxes)");
        }

        eEqnViscous.relax();

        fvOptions.constrain(eEqnViscous);

        eEqnViscous.solve();

        fvOptions.correct(e);

        et.ref() = e() - ev();

        e.correctBoundaryConditions();
        et.boundaryFieldRef() = e.boundaryField() - ev.boundaryField();
    }
    else
    {
        // Two-temperature models: trans-rotational operator splitting —
        // et diffused with kappatr, ev handled by the ev-equation.
        // See README: fixes the ~19% translational wall heat flux error.
        fvScalarMatrix eEqnViscous
        (
            fvm::ddt(rho, et) - fvc::ddt(rho, et)
          - fvm::laplacian(thermo.kappatr()/thermo.CvtrMix(), et)
          ==
            fvOptions(rho, et)
        );

        if (speciesDiffusionName != "noSpeciesDiffusion")
        {
            if (downgradeSingleTv)
            {
                // Two-temperature model, single Tv: the et-equation must only
                // receive the trans-rotational share of the species diffusion
                // enthalpy flux, sum_j hts_j*J_j. The vibro-electronic share,
                // returned by multiSpeciesVEHeatSource(), is already transported
                // by the ev-equation (evEqnViscous_SingleTv.H), so it is
                // subtracted here. Identity used: hs = hts + hevel, with
                // hvs == evs and hels == eels (the pv share of the enthalpy
                // belongs to the translational mode only).
                eEqnViscous +=
                    fvc::div
                    (
                        speciesDiffusion().multiSpeciesHeatSource()
                      - speciesDiffusion().multiSpeciesVEHeatSource(),
                        "div(diffusiveFluxes)"
                    );

                // Catalytic wall: boundary completion of the
                // trans-rotational share
                // 催化壁面：平动-转动份额的边界完成项
                if (catalyticWallPresent)
                {
                    eEqnViscous +=
                        fvc::div
                        (
                            catalyticFluxTr(),
                            "div(diffusiveFluxes)"
                        );
                }
            }
            else
            {
                // Two-temperature model, one Tv per molecule: mirror exactly the
                // vibro-electronic share that the per-species ev-equations
                // (evEqnViscous.H) receive, sum_k hevel_k*J_k over the
                // solved-vibrational-equation species, and subtract it from the
                // full species diffusion enthalpy flux so that the et-equation
                // is left with the trans-rotational share only.
                volVectorField veShare
                (
                    IOobject
                    (
                        "veShare",
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedVector
                    (
                        "zero",
                        dimEnergy/dimArea/dimTime,
                        vector::zero
                    )
                );

                for (label speciek = 0; speciek < noSolvedVibEq; speciek++)
                {
                    veShare += speciesDiffusion().Jcorr(speciek)*hevk[speciek];
                }

                eEqnViscous +=
                    fvc::div
                    (
                        speciesDiffusion().multiSpeciesHeatSource() - veShare,
                        "div(diffusiveFluxes)"
                    );

                // Catalytic wall: boundary completion of the
                // trans-rotational share
                // 催化壁面：平动-转动份额的边界完成项
                if (catalyticWallPresent)
                {
                    eEqnViscous +=
                        fvc::div
                        (
                            catalyticFluxTr(),
                            "div(diffusiveFluxes)"
                        );
                }
            }
        }

        eEqnViscous.relax();

        fvOptions.constrain(eEqnViscous);

        eEqnViscous.solve();

        fvOptions.correct(et);

        e.ref() = et() + ev();

        et.correctBoundaryConditions();
        e.boundaryFieldRef() = et.boundaryField() + ev.boundaryField();
    }

    rhoE = rho*(e + 0.5*magSqr(U));
}

if (hyLight)
{
    thermo.correct2T_Light();
}
else
{
    thermo.correct2T();
    thermo.correctOverallTemperature();
}
