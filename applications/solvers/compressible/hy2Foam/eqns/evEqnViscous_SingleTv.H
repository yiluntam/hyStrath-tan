
// --- Solve vibrational energy equation, VISCOUS
if (!inviscid)
{
    fvScalarMatrix evEqnViscous
    (
        fvm::ddt(rho, ev) - fvc::ddt(rho, ev)
      - fvm::laplacian(transport().kappave()/thermo.CvvelMix(), ev)
      ==
        fvOptions(rho, ev)
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        evEqnViscous +=
            fvc::div
            (
                speciesDiffusion().multiSpeciesVEHeatSource(),
                "div(diffusiveFluxes)"
            );
    }

    // Catalytic wall: boundary completion of the vibro-electronic share
    // (the interior source term has zero boundary values)
    // 催化壁面：振电份额的边界完成项
    if (catalyticWallPresent)
    {
        evEqnViscous +=
            fvc::div
            (
                catalyticFluxVe(),
                "div(diffusiveFluxes)"
            );
    }

    evEqnViscous.solve();

    ev.correctBoundaryConditions();
    rhoEv = rho*ev;
}

