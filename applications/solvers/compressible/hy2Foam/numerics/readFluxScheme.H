
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
fluxMap["AUSM"] = AUSMFlux;
fluxMap["Roe"] = RoeFlux;

// --- 两级配置 ---
//   fluxScheme : 通量格式  = Kurganov (原生 KNP, 默认) | AUSM | Roe
//   convSchemes: 对流重构  = MUSCL | WENOZ5 | WENOZ7 (仅 AUSM/Roe 生效)
word fluxScheme("Kurganov");
word convSchemes("MUSCL");
FixedList<word, 2> blendedSchemes;

if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    // 旧组合名自动映射到新的两级语法 (仅提示)
    if (fluxScheme == "MUSCL_AUSM")
    {
        fluxScheme = "AUSM";
        convSchemes = "MUSCL";
        Info<< "fluxScheme MUSCL_AUSM is deprecated -> "
            << "fluxScheme AUSM + convSchemes MUSCL" << endl;
    }
    else if (fluxScheme == "WENOZ_AUSM")
    {
        fluxScheme = "AUSM";
        convSchemes = "WENOZ5";
        Info<< "fluxScheme WENOZ_AUSM is deprecated -> "
            << "fluxScheme AUSM + convSchemes WENOZ5" << endl;
    }
    else if (fluxScheme == "WENOZ7_AUSM")
    {
        fluxScheme = "AUSM";
        convSchemes = "WENOZ7";
        Info<< "fluxScheme WENOZ7_AUSM is deprecated -> "
            << "fluxScheme AUSM + convSchemes WENOZ7" << endl;
    }

    if
    (
        (fluxScheme == "Kurganov")
     or (fluxScheme == "AUSM")
     or (fluxScheme == "Roe")
     or (fluxScheme == "Tadmor")
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
            << " is not a valid choice. Options are: Kurganov, AUSM, "
            << "Roe, Tadmor, and blended. The reconstruction is "
            << "selected separately via 'convSchemes' (MUSCL, WENOZ5, "
            << "WENOZ7) together with fluxScheme AUSM or Roe."
            << abort(FatalError);
    }
}

// --- AUSM+up 远场马赫数 Mach_inf (低马赫修正参数 fa/Mref 需要 M∞;
//     仅 fluxScheme = AUSM 时必须为正, 其余格式指定了也不起作用) ---
if (fluxScheme == "AUSM")
{
    const scalar MachInf =
        mesh.schemesDict().lookupOrDefault<scalar>("Mach_inf", -1.0);

    if (MachInf <= 0.0)
    {
        FatalErrorIn("readFluxScheme")
            << "AUSM+up requires a positive freestream Mach number: "
            << "please set 'Mach_inf' (the case freestream Mach number) "
            << "in fvSchemes."
            << abort(FatalError);
    }

    ausmMachInfSq = sqr(MachInf);
    Info<< "AUSM+up Mach_inf = " << MachInf << endl;
}
else if (mesh.schemesDict().found("Mach_inf"))
{
    Info<< "Mach_inf specified but ignored for fluxScheme "
        << fluxScheme << endl;
}

// --- convSchemes: 对流重构格式 (仅 fluxScheme = AUSM/Roe 时读取并生效) ---
// 原生 KNP/Tadmor 路径不读取任何新选项, 数值结果与原 hyStrath 完全一致
if (fluxScheme == "AUSM" or fluxScheme == "Roe")
{
    mesh.schemesDict().readIfPresent("convSchemes", convSchemes);

    if
    (
        convSchemes != "MUSCL"
     and convSchemes != "WENOZ5"
     and convSchemes != "WENOZ7"
    )
    {
        FatalErrorIn("readFluxScheme")
            << "convSchemes: " << convSchemes
            << " is not a valid choice. Options are: MUSCL, WENOZ5, "
            << "and WENOZ7."
            << abort(FatalError);
    }
    Info<< "convSchemes: " << convSchemes << endl;

    // --- MUSCL 限制器与 kappa 系数 (仅 convSchemes = MUSCL 时读取并生效;
    //     WENOZ5/WENOZ7 不使用, 指定了也不起作用) ---
    if (convSchemes == "MUSCL")
    {
        word musclLimiter("vanAlbada");
        if
        (
            mesh.schemesDict().readIfPresent("MUSCLLimiter", musclLimiter)
         and musclLimiter != "minmod"
         and musclLimiter != "vanLeer"
         and musclLimiter != "superbee"
         and musclLimiter != "vanAlbada"
        )
        {
            FatalErrorIn("readFluxScheme")
                << "MUSCLLimiter: " << musclLimiter
                << " is not a valid choice. Options are: minmod, vanLeer, "
                << "superbee, and vanAlbada."
                << abort(FatalError);
        }
        musclLimiterType = (musclLimiter == "minmod") ? 0 :
                           (musclLimiter == "vanLeer") ? 1 :
                           (musclLimiter == "superbee") ? 2 : 3;
        musclKappa =
            mesh.schemesDict().lookupOrDefault<scalar>("MUSCLKappa", 1.0/3.0);

        Info<< "MUSCL limiter: " << musclLimiter
            << ", kappa = " << musclKappa << endl;
    }
}
else if (mesh.schemesDict().found("convSchemes"))
{
    Info<< "convSchemes specified but ignored for fluxScheme "
        << fluxScheme << endl;
}

// --- Roe 熵修正系数 EntropyCoeff (Harten-Hyman 1983 逐波熵修正,
//     仅 fluxScheme = Roe 时读取并生效, 其他格式指定了也不起作用) ---
//   ε_k = EntropyCoeff*4*max(0, max(λ̂_k - λ_k,L, λ_k,R - λ̂_k)),
//   λ_k,L/R 为该波在面两侧状态的特征值 (SU2 NEMO roe.cpp);
//   |λ̂| < ε: (λ̂²+ε²)/(2ε), 否则 |λ̂|。默认 1.0 = SU2 原样;
//   波不跨越 (ε = 0) 时修正自动关闭。
if (fluxScheme == "Roe")
{
    roeEntropyFixCoeff =
        mesh.schemesDict().lookupOrDefault<scalar>("EntropyCoeff", 1.0);

    Info<< "Roe Harten-Hyman entropy fix: EntropyCoeff = "
        << roeEntropyFixCoeff << endl;
}
else if (mesh.schemesDict().found("EntropyCoeff"))
{
    Info<< "EntropyCoeff specified but ignored for fluxScheme "
        << fluxScheme << endl;
}


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

// --- 处理器边界 ghost 层 (并行 MUSCL/WENO 模板用) ---
// ghostStencilBuild.H 初始化, ghostFieldExchange.H (经
// MUSCL_WENO_Reconstruct.H) 每时间步使用; 串行时 nGhostCells = 0
label ghostDepth = 0;           // ghost 层数: MUSCL=1, WENOZ5=2, WENOZ7=3
label nGhostCells = 0;          // ghost 单元总数 (串行 = 0)
labelList procPatchIDs;         // 参与交换的 processor patch
labelListList ghostSendCells;   // 每 patch 发送单元表 [face*ghostDepth + d-1]
labelList ghostStart;           // 每 patch 的 ghost 编号起始偏移

// 仅当选用 AUSM/Roe (需 convSchemes 重构) 时初始化; KNP/KT 不需要
if
(
    fluxScheme == "AUSM"
 or fluxScheme == "Roe"
)
{
    #include "structuredStencilBuild.H"
    #include "ghostStencilBuild.H"
}

