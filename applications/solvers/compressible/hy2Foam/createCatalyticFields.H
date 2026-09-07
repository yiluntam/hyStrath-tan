//-----------------------------------------------------------
// Catalytic wall boundary conditions: patch detection,
// start-up validation and completion-flux fields
// 催化壁面边界条件：壁面检测、启动校验与能量完成通量场
//-----------------------------------------------------------

// Detection is done by runtime type name only (pure word comparison),
// so the solver has no link-time dependency on the BC library:
// 壁面检测仅按运行时类型名（纯字符串比较），求解器不引入对
// BC 库的链接期依赖：
//   nonEqScottCatalyticY        finite-rate recombination, atoms 原子
//   nonEqScottCatalyticProductY recombination product 产物
//   nonEqSuperCatalyticY        supercatalytic, prescribed Y_w 超催化
// Without such BCs, catalyticWallPresent stays false and every
// solver-side addition is skipped (regression invariant).
// 无此类边界条件时 catalyticWallPresent 保持 false，求解器侧的全部
// 改动被跳过（回归不变量）。

bool catalyticWallPresent = false;
labelList catalyticPatchIDs;

// Completion flux fields for the energy equations, constructed only when
// catalytic walls are present; their boundary values are refreshed every
// iteration in eqns/catalyticWallSources.H (included at the end of
// YEqn.H). All fields are [W/m^2] vector fields:
// 能量方程的完成通量场，仅在存在催化壁面时构造；其边界面值在
// YEqn.H 末尾的 eqns/catalyticWallSources.H 中逐迭代刷新：
//   catalyticFluxE   single-T mode, full sensible share   -> e-equation
//   catalyticFluxTr  two-T modes, trans-rotational share  -> et-equation
//   catalyticFluxVe  single-Tv mode, vibro-electronic     -> ev-equation
//   catalyticFluxVeK per-molecule-Tv, per-species share   -> evk-equations
autoPtr<volVectorField> catalyticFluxE;
autoPtr<volVectorField> catalyticFluxTr;
autoPtr<volVectorField> catalyticFluxVe;
autoPtr<PtrList<volVectorField> > catalyticFluxVeK;

// Reported catalytic wall heat flux, absolute-enthalpy basis, positive
// into the wall (same convention as the other wallHeatFlux fields).
// Written explicitly from write.H (NO_WRITE here).
// 报告用催化壁面热流（绝对焓基准，指向壁面为正，与其它壁面热流
// 场同约定）；由 write.H 显式写出。
volScalarField wallHeatFlux_cat
(
    IOobject
    (
        "wallHeatFlux_cat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("wallHeatFlux_cat", dimensionSet(1, 0, -3, 0, 0), 0.0)
);

{
    DynamicList<label> catalyticPatchIDsList;

    forAll(Y, speciei)
    {
        forAll(Y[speciei].boundaryField(), patchi)
        {
            const word& bType =
                Y[speciei].boundaryField()[patchi].type();

            if
            (
                bType == "nonEqScottCatalyticY"
             || bType == "nonEqScottCatalyticProductY"
             || bType == "nonEqSuperCatalyticY"
            )
            {
                if (findIndex(catalyticPatchIDsList, patchi) == -1)
                {
                    catalyticPatchIDsList.append(patchi);
                }
            }
        }
    }

    catalyticPatchIDs = catalyticPatchIDsList;
    catalyticWallPresent = (catalyticPatchIDs.size() > 0);
}

if (catalyticWallPresent)
{
    Info<< "\nCatalytic wall patches detected:"
        << catalyticPatchIDs << endl;

    if (speciesDiffusionName == "noSpeciesDiffusion")
    {
        FatalErrorIn("createCatalyticFields.H")
            << "Catalytic wall boundary conditions require a multi-species"
            << " diffusion model. Set 'multiSpeciesTransport' (e.g. to"
            << " Fick) in the 'transportModels' sub-dictionary of"
            << " constant/transportProperties."
            << exit(FatalError);
    }

    if (!solveSpeciesEqns)
    {
        FatalErrorIn("createCatalyticFields.H")
            << "Catalytic wall boundary conditions require the species"
            << " transport equations to be solved."
            << exit(FatalError);
    }

    if (inviscid)
    {
        WarningIn("createCatalyticFields.H")
            << "Catalytic wall boundary conditions have no effect on the"
            << " energy budget in inviscid mode: the species-diffusive"
            << " enthalpy flux is not part of the inviscid energy"
            << " equation. The species wall conditions themselves remain"
            << " active." << endl;
    }

    // Per-patch consistency checks 逐壁面一致性检查
    forAll(catalyticPatchIDs, i)
    {
        const label patchi = catalyticPatchIDs[i];
        const word& patchName = mesh.boundaryMesh()[patchi].name();

        bool hasScott = false;
        bool hasProduct = false;
        bool hasSuper = false;
        label nSuper = 0;
        scalar superSum = 0.0;

        forAll(Y, speciei)
        {
            const word& bType =
                Y[speciei].boundaryField()[patchi].type();

            if (bType == "nonEqScottCatalyticY")
            {
                hasScott = true;
            }
            else if (bType == "nonEqScottCatalyticProductY")
            {
                hasProduct = true;
            }
            else if (bType == "nonEqSuperCatalyticY")
            {
                hasSuper = true;
                nSuper++;
                superSum += gSum(Y[speciei].boundaryField()[patchi]);
            }
        }

        if (hasScott && hasSuper)
        {
            FatalErrorIn("createCatalyticFields.H")
                << "Patch '" << patchName << "': finite-rate catalytic"
                << " (nonEqScottCatalyticY) and supercatalytic"
                << " (nonEqSuperCatalyticY) boundary conditions are"
                << " mutually exclusive on the same patch."
                << exit(FatalError);
        }

        if (hasScott && !hasProduct)
        {
            WarningIn("createCatalyticFields.H")
                << "Patch '" << patchName << "': nonEqScottCatalyticY is"
                << " set on (an) atomic species but no species uses"
                << " nonEqScottCatalyticProductY on this patch; the"
                << " recombined mass cannot re-enter the gas phase and"
                << " the wall mass balance is not closed." << endl;
        }

        if (hasSuper)
        {
            Info<< "Catalytic wall: patch '" << patchName
                << "' is supercatalytic; sum of the prescribed wall mass"
                << " fractions = " << superSum << " (" << nSuper
                << " species)" << endl;

            if (nSuper == Y.size() && mag(superSum - 1.0) > 1.0e-6)
            {
                WarningIn("createCatalyticFields.H")
                    << "Patch '" << patchName << "': all species are"
                    << " prescribed (supercatalytic) and the wall mass"
                    << " fractions sum to " << superSum << ", not 1."
                    << endl;
            }
        }
    }

    // Construct the completion flux fields 构造完成通量场
    const dimensionSet dimCatFlux(dimEnergy/dimArea/dimTime);

    catalyticFluxE.set
    (
        new volVectorField
        (
            IOobject
            (
                "catalyticFluxE",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimCatFlux, vector::zero)
        )
    );

    catalyticFluxTr.set
    (
        new volVectorField
        (
            IOobject
            (
                "catalyticFluxTr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimCatFlux, vector::zero)
        )
    );

    catalyticFluxVe.set
    (
        new volVectorField
        (
            IOobject
            (
                "catalyticFluxVe",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("zero", dimCatFlux, vector::zero)
        )
    );

    catalyticFluxVeK.set(new PtrList<volVectorField>(noSolvedVibEq));

    for (label speciek = 0; speciek < noSolvedVibEq; speciek++)
    {
        catalyticFluxVeK->set
        (
            speciek,
            new volVectorField
            (
                IOobject
                (
                    "catalyticFluxVeK_"
                  + composition.solvedVibEqSpecies()[speciek],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimCatFlux, vector::zero)
            )
        );
    }
}
