
// --- Solve vibrational energy equations, VISCOUS
if (!inviscid)
{
    forAll(rhoEvk, speciei)
    {
        DynamicList<label> removedCellsFromMatrix;
        DynamicList<scalar> valuesToImpose;

        forAll(Y[speciei], celli)
        {
            if
            (
                Y[speciei][celli]
                    < Foam::rho2ReactionThermo::miniYforSolvingEvEqn
            )
            {
                removedCellsFromMatrix.append(celli);
            }
        }

        forAll(removedCellsFromMatrix, celli)
        {
            valuesToImpose.append(0.0);
        }

        fvScalarMatrix evEqnViscous
        (
            fvm::ddt(pD[speciei], hevk[speciei])
          - fvc::ddt(pD[speciei], hevk[speciei])
          - fvm::laplacian(transport().alphave(speciei), hevk[speciei])
        );

        if (speciesDiffusionName != "noSpeciesDiffusion")
        {
            evEqnViscous +=
                fvc::div
                (
                    speciesDiffusion().Jcorr(speciei)*hevk[speciei],
                    "div(diffusiveFluxes)"
                );
        }

        // Catalytic wall: boundary completion of the vibro-electronic
        // share of this species (the interior source term has zero
        // boundary values)
        // 催化壁面：本组分振电份额的边界完成项
        if (catalyticWallPresent)
        {
            evEqnViscous +=
                fvc::div
                (
                    catalyticFluxVeK()[speciei],
                    "div(diffusiveFluxes)"
                );
        }

        evEqnViscous.setValues(removedCellsFromMatrix, valuesToImpose);

        evEqnViscous.solve();

        rhoEvk[speciei] = pD[speciei]*hevk[speciei];
    }

    thermo.correctHEVK();

    rhoEv *= 0.0;
    forAll(hevk, speciei)
    {
        rhoEv += pD[speciei]*hevk[speciei];
    }

    ev = rhoEv/rho;
}

