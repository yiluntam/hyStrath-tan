
// --- Solve translational-rotational energy, VISCOUS
if (!inviscid)
{
    fvScalarMatrix eEqnViscous
    (
        fvm::ddt(rho, et) - fvc::ddt(rho, et)
      - fvm::laplacian(thermo.kappatr()/thermo.CvtrMix(), et)
      ==
        fvOptions(rho, et)
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        if (downgradeSingleT)
        {
            // Single-temperature model: no separate ev-equation is solved,
            // so the full species diffusion enthalpy flux sum_j hs_j*J_j
            // enters the energy budget here.
            eEqnViscous +=
                fvc::div
                (
                    speciesDiffusion().multiSpeciesHeatSource(),
                    "div(diffusiveFluxes)"
                );
        }
        else if (downgradeSingleTv)
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
        }
    }

    eEqnViscous.relax();

    fvOptions.constrain(eEqnViscous);

    eEqnViscous.solve();

    fvOptions.correct(et);

    e.ref() = et() + ev();

    if (downgradeSingleT)
    {
        e.correctBoundaryConditions();
        et.boundaryFieldRef() = e.boundaryField() - ev.boundaryField();
    }
    else
    {
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
