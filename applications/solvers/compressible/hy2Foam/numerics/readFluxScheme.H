
typedef void (*FnPtr)
(
    const bool&,
    const bool&,
    const bool&,
    const scalar&,
    const vector&,
    const scalar&,
    const scalar&,
    const scalar&,
    const scalar&,
    const PtrList<scalar>&,
    const PtrList<scalar>&,
    const vector&,
    const vector&,
    const vector&,
    const vector&,
    const scalar&,
    const scalar&,
    const PtrList<scalar>&,
    const PtrList<scalar>&,
    const scalar&,
    const scalar&,
    scalar&,
    vector&,
    scalar&,
    vector&,
    PtrList<scalar>&,
    scalar&
);

std::map<word, FnPtr> fluxMap;
fluxMap["Kurganov"] = KurganovFlux;
fluxMap["Tadmor"] = TadmorFlux;
fluxMap["MUSCL_AUSM"] = AUSMFlux;
fluxMap["WENOZ_AUSM"] = AUSMFlux;
fluxMap["WENOZ7_AUSM"] = AUSMFlux;

word fluxScheme("Kurganov");
FixedList<word, 2> blendedSchemes;

if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    if
    (
        (fluxScheme == "Tadmor")
     or (fluxScheme == "Kurganov")
     or (fluxScheme == "MUSCL_AUSM")
     or (fluxScheme == "WENOZ_AUSM")
     or (fluxScheme == "WENOZ7_AUSM")
    )
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else if (fluxScheme == "blended")
    {
        if (mesh.schemesDict().readIfPresent("blendedSchemes", blendedSchemes))
        {
            Info<< "fluxScheme for shock vicinity: " << blendedSchemes[0]
                << endl;
            Info<< "fluxScheme for near-wall region: " << blendedSchemes[1]
                << endl;
        }
        else
        {
            FatalErrorIn
            (
                "readFluxScheme"
            )   << "blendedSchemes: " << fluxScheme
                << " is not a valid choice. "
                << "Options are: Kurganov, and Tadmor."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "readFluxScheme"
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. Options are: Kurganov, Tadmor, "
            << "MUSCL_AUSM, WENOZ_AUSM, WENOZ7_AUSM, and blended"
            << abort(FatalError);
    }
}

scalar sqrMachInf =
    sqr(mesh.schemesDict().lookupOrDefault<scalar>("Mach_inf", SMALL));


// --- Declarations ---
surfaceScalarField amaxSf
(
    IOobject
    (
        "amaxSf",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("amaxSf", dimVolume/dimTime, 0.0)
);

surfaceVectorField aByU
(
    IOobject
    (
        "aByU",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("aByU", dimVelocity, vector::zero)
);

surfaceVectorField phiUp
(
    IOobject
    (
        "phiUp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("phiUp", dimMass/dimTime*dimVelocity, vector::zero)
);

surfaceScalarField phiEp
(
    IOobject
    (
        "phiEp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("phiEp", dimEnergy/dimTime, 0.0)
);

PtrList<surfaceScalarField> phiEvk(noSolvedVibEq);

if (downgradeSingleT)
{
    phiEvk.clear();
}
else if (downgradeSingleTv)
{
    phiEvk.resize(1);

    phiEvk.set
    (
        0,
        new surfaceScalarField
        (
            IOobject
            (
                "phiEvk",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("phiEvk", dimEnergy/dimTime, 0.0)
        )
    );
}
else
{
    forAll(phiEvk, speciei)
    {
        phiEvk.set
        (
            speciei,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiEvk" + Y[speciei].name(),
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "phiEvk" + Y[speciei].name(),
                    dimEnergy/dimTime,
                    0.0
                )
            )
        );
    }
}

// --- 结构化网格模板 (MUSCL/WENO 高精度重构用) ---
// 变量在 if 外部声明, 使其在函数作用域内对后续 include (MUSCL_WENO_Reconstruct.H) 可见
labelList axisFaces;
label nAxisFaces = 0;
labelList cell_plus_x, cell_minus_x;
labelList cell_plus_y, cell_minus_y;
labelList cell_plus_z, cell_minus_z;
bool isStructuredGrid = false;

// 仅当选用 MUSCL/WENO 格式时初始化; KNP/KT 不需要
if
(
    fluxScheme == "MUSCL_AUSM"
 or fluxScheme == "WENOZ_AUSM"
 or fluxScheme == "WENOZ7_AUSM"
)
{
    #include "structuredStencilBuild.H"
}

