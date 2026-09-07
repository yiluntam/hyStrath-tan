//------------------------------------------------------------------------------
// AUSM+-up flux function (Liou, 2006) — SU2 CUpwAUSMPLUSUP_Flow 对齐
//
// - 界面声速: 临界声速 min 形式 (逐侧 γ, NEMO AUSM+up2 同法)
//     a*L = √(2(γL-1)/(γL+1)·HL),  âL = a*L²/max(a*L, u_nL)
//     a*R = √(2(γR-1)/(γR+1)·HR),  âR = a*R²/max(a*R, -u_nR)
//     aF = min(âL, âR)
//   γ 由 c² = γp/ρ 反解 (两温逐侧, 与 Roe_Flux.H 一致)。
// - 低马赫参数 (M∞ 由 fvSchemes 的 Mach_inf 指定, AUSM 时必须为正):
//     MF² = ½(mL² + mR²),  Mref² = min(1, max(MF², M∞²))
//     fa = 2√Mref² - Mref²,  α = 3/16*(-4 + 5fa²)
// - M±(4) (β = 1/8) 与压力分裂系数 β± = P±(5) (动态 α):
//     Mp = -(Kp/fa)*max(1-MF², 0)*Δp/(ρ̄*aF²),  Kp = 0.25
//     Pu = -Ku*fa*βL*βR*2ρ̄*aF*Δu_n,            Ku = 0.75
//   质量通量 ṁ = aF*(max(mF,0)ρL + min(mF,0)ρR), mF = M+L + M-R + Mp;
//   界面压力 P½ = βL*pL + βR*pR + Pu。
//   M∞ ≥ 1 时 fa = 1, α = 3/16, 退化为 Liou 2006 高马赫极限。
// - 振动能量按被动标量随迎风质量通量输运 (ṁ·Y·ev);
//   ALE: u_n 相对网格速度。
// - 函数签名与 KurganovFlux / TadmorFlux 一致, 注册于同一 fluxMap。
//------------------------------------------------------------------------------

// --- AUSM+up 远场马赫数平方 (由 readFluxScheme.H 从 fvSchemes 的
//     Mach_inf 设置; 默认 1 对应 M∞ ≥ 1 的 fa = 1 高马赫极限) ---
double ausmMachInfSq = 1.0;

void AUSMFlux
(
    const bool& downgradeSingleT,
    const bool& downgradeSingleTv,
    const bool& meshMoving,
    const scalar& mesh_phi,
    const vector& Sf,
    const scalar& cSf_pos,
    const scalar& cSf_neg,
    const scalar& rho_pos,
    const scalar& rho_neg,
    const PtrList<scalar>& pD_pos,
    const PtrList<scalar>& pD_neg,
    const vector& U_pos,
    const vector& U_neg,
    const vector& rhoU_pos,
    const vector& rhoU_neg,
    const scalar& p_pos,
    const scalar& p_neg,
    const PtrList<scalar>& evk_pos,
    const PtrList<scalar>& evk_neg,
    const scalar& e_pos,
    const scalar& e_neg,
    scalar& amaxSf,
    vector& aByU,
    scalar& rhoFlux,
    vector& rhoUFlux,
    PtrList<scalar>& rhoEvkFlux,
    scalar& rhoEFlux
)
{
    // --- 面单位法向与面积 ---
    const scalar magSf = mag(Sf);
    const vector n = Sf / magSf;

    // --- 面投影速度 (ALE: 相对网格) ---
    scalar un_pos = U_pos & n;
    scalar un_neg = U_neg & n;

    if (meshMoving)
    {
        const scalar mesh_un = mesh_phi / magSf;
        un_pos -= mesh_un;
        un_neg -= mesh_un;
    }

    // --- 每侧冻结比热比 (c² = γp/ρ 反解; 两温逐侧, 与 Roe_Flux.H 一致) ---
    const scalar c_pos = cSf_pos / magSf;
    const scalar c_neg = cSf_neg / magSf;
    const scalar gammaL = rho_pos*sqr(c_pos)/max(p_pos, SMALL);
    const scalar gammaR = rho_neg*sqr(c_neg)/max(p_neg, SMALL);

    // --- 界面声速 (SU2: 临界声速 min 形式) ---
    // H = e + p/ρ + ½|U|² (e 为比内能, hy2Foam 约定)
    const scalar H_pos = e_pos + p_pos/rho_pos + 0.5*magSqr(U_pos);
    const scalar H_neg = e_neg + p_neg/rho_neg + 0.5*magSqr(U_neg);

    const scalar aStarL =
        Foam::sqrt(max(2.0*(gammaL - 1.0)/(gammaL + 1.0)*H_pos, 0.0));
    const scalar aStarR =
        Foam::sqrt(max(2.0*(gammaR - 1.0)/(gammaR + 1.0)*H_neg, 0.0));

    const scalar aHatL = sqr(aStarL)/max(aStarL, un_pos);
    const scalar aHatR = sqr(aStarR)/max(aStarR, -un_neg);
    const scalar aF = max(min(aHatL, aHatR), SMALL);

    // --- 马赫数与低马赫参数 (SU2: MF², Mref², fa, α) ---
    const scalar mL = un_pos / aF;
    const scalar mR = un_neg / aF;

    const scalar MFsq = 0.5*(sqr(mL) + sqr(mR));
    const scalar Mrefsq = min(1.0, max(MFsq, ausmMachInfSq));
    const scalar fa = 2.0*Foam::sqrt(Mrefsq) - Mrefsq;
    const scalar alpha = 3.0/16.0*(-4.0 + 5.0*fa*fa);
    const scalar beta = 1.0/8.0;

    // --- M±(4) 分裂与压力分裂系数 β± = P±(5) (动态 α) ---
    scalar mLP, betaLP, mRM, betaRM;

    if (mag(mL) <= 1.0)
    {
        const scalar p1 = 0.25*sqr(mL + 1.0);
        const scalar p2 = sqr(sqr(mL) - 1.0);
        mLP = p1 + beta*p2;
        betaLP = p1*(2.0 - mL) + alpha*mL*p2;
    }
    else
    {
        mLP = 0.5*(mL + mag(mL));
        betaLP = mLP/mL;
    }

    if (mag(mR) <= 1.0)
    {
        const scalar p1 = 0.25*sqr(mR - 1.0);
        const scalar p2 = sqr(sqr(mR) - 1.0);
        mRM = -p1 - beta*p2;
        betaRM = p1*(2.0 + mR) - alpha*mR*p2;
    }
    else
    {
        mRM = 0.5*(mR - mag(mR));
        betaRM = mRM/mR;
    }

    // --- 压力/速度扩散项 (Kp = 0.25, σ = 1, Ku = 0.75) ---
    const scalar Kp = 0.25;
    const scalar sigma = 1.0;
    const scalar Ku = 0.75;
    const scalar rhoF = 0.5*(rho_pos + rho_neg);

    const scalar MpDiss = -(Kp/max(fa, SMALL))
        *max(1.0 - sigma*MFsq, 0.0)
        *(p_neg - p_pos)/(rhoF*sqr(aF) + SMALL);

    const scalar Pu = -Ku*fa*betaLP*betaRM*2.0*rhoF*aF*(un_neg - un_pos);

    // --- 界面马赫数与压力 ---
    const scalar mF = mLP + mRM + MpDiss;
    const scalar P_half = betaLP*p_pos + betaRM*p_neg + Pu;

    // --- 质量通量 (SU2: ṁ = aF*(max(mF,0)ρL + min(mF,0)ρR)) ---
    const scalar mDot = aF*(max(mF, 0.0)*rho_pos + min(mF, 0.0)*rho_neg);

    // --- 迎风状态 (由 mF 符号选择) ---
    vector U_upwind;   // upwind velocity
    scalar e_upwind;   // upwind total specific energy
    scalar p_upwind;   // upwind pressure (for enthalpy)
    scalar rho_upwind; // upwind density

    if (mF > 0.0)
    {
        // Flow from left (pos) to right (neg)
        rho_upwind = rho_pos;
        U_upwind   = U_pos;
        e_upwind   = e_pos;
        p_upwind   = p_pos;
    }
    else
    {
        // Flow from right (neg) to left (pos)
        rho_upwind = rho_neg;
        U_upwind   = U_neg;
        e_upwind   = e_neg;
        p_upwind   = p_neg;
    }

    // Total enthalpy H = e + p/rho + 0.5*|U|^2
    // (e is the INTERNAL energy in hy2Foam, cf. eEqnInviscid.H:23
    //  e = rhoE/rho - 0.5*magSqr(U); KNP.H adds the same kinetic term)
    const scalar H_upwind = e_upwind + p_upwind/rho_upwind
                          + 0.5*magSqr(U_upwind);

    // --- Face fluxes ---

    // Density flux: ṁ * |Sf|
    rhoFlux = mDot * magSf;

    // Momentum flux: ṁ * U_upwind * |Sf| + P_half * Sf
    rhoUFlux = mDot * U_upwind * magSf + P_half * Sf;

    // Energy flux: ṁ * H_upwind * |Sf|
    // In mesh-moving case, add mesh_phi * P_half for pressure work
    rhoEFlux = mDot * H_upwind * magSf;
    if (meshMoving)
    {
        rhoEFlux += mesh_phi * P_half;
    }

    // Vibrational energy fluxes
    if (downgradeSingleT)
    {
        rhoEvkFlux.clear();
    }
    else if (downgradeSingleTv)
    {
        // Single vib. mode: ṁ * ev_upwind * |Sf|
        if (mF > 0.0)
        {
            rhoEvkFlux[0] = mDot * evk_pos[0] * magSf;
        }
        else
        {
            rhoEvkFlux[0] = mDot * evk_neg[0] * magSf;
        }
    }
    else
    {
        // Per-species vib. energy: ṁ * Y_s * ev_s * |Sf|
        // ṁ = aF*(...)*rho_upwind, Y_s = pD_s/rho, so the rho factors cancel:
        // ṁ * Y * evk * |Sf| = aF*(...)*pD * evk * |Sf|  [J/s]
        forAll(rhoEvkFlux, speciei)
        {
            if (mF > 0.0)
            {
                rhoEvkFlux[speciei] = mDot * (pD_pos[speciei]/rho_pos) * evk_pos[speciei] * magSf;
            }
            else
            {
                rhoEvkFlux[speciei] = mDot * (pD_neg[speciei]/rho_neg) * evk_neg[speciei] * magSf;
            }
        }
    }

    // --- amaxSf for CFL: (aF + max(|un_L|, |un_R|)) * |Sf| ---
    amaxSf = (aF + max(mag(un_pos), mag(un_neg))) * magSf;

    // --- aByU for viscous stress (σ·U term) ---
    aByU = U_upwind;
}
