
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
        eEqnViscous +=
            fvc::div
            (
                speciesDiffusion().multiSpeciesHeatSource(),
                "div(diffusiveFluxes)"
            );
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
