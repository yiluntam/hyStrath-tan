//------------------------------------------------------------------------------
// Roe 通量函数 (Roe, 1981) + Harten-Hyman 熵修正 (Harten & Hyman, 1983)
//
// - Roe 平均: ρ̂ = √(ρL·ρR), Û 与 Ĥ 为 √ρ 加权平均
//   (H = e + p/ρ + ½|U|², e 为比内能, hy2Foam 约定, 参见 eEqnInviscid.H)
// - 声速: ĉ² = (γ̂-1)(Ĥ - ½|Û|² - ê_v)
//   γ̂ 为两侧冻结声速反解值 (γ = ρc²/p) 的 √ρ 加权平均, ê_v 为 √ρ 加权
//   平均比振动能。单温完全气体 (ê_v = 0, γ 为常数) 时退化为经典
//   ĉ² = (γ̂-1)(Ĥ - ½|Û|²)。两温气体不存在严格 Roe 线性化, 此为与
//   求解器冻结平动声速 (c = √(gammatr·rPsi)) 一致的近似。
// - 耗散项 D = Σ_k |λ̂_k| α_k K_k, 波强
//     α₁ = (Δp - ρ̂ĉΔu_n)/(2ĉ²)   (左声波)
//     α₂ = Δρ - Δp/ĉ²              (熵/接触波)
//     α₃ = (Δp + ρ̂ĉΔu_n)/(2ĉ²)   (右声波)
//     β_i = ρ̂Δu_ti                 (剪切波, 并入接触特征值)
//   闭式代数已与解析 Jacobian 特征分解逐项核对 (残差 ~1e-12)。
// - 熵修正 (Harten-Hyman 1983 分段二次型, 逐波带宽, SU2 NEMO roe.cpp):
//     |λ̂| ≥ ε: |λ̂| 不变
//     |λ̂| < ε: (λ̂²+ε²)/(2ε)   (C1 连续, 恒有 ≥ |λ̂|)
//   ε_k = EntropyCoeff*4*max(0, max(λ̂_k - λ_k,L, λ_k,R - λ̂_k)) 由
//   该波两侧状态特征值张成 (λ_k,L/R 为面两侧状态的同类波特征值);
//   系数由 fvSchemes 的 EntropyCoeff 指定 (默认 1.0 = SU2 原样,
//   仅 fluxScheme = Roe 时生效, 见 readFluxScheme.H)。波不跨越时
//   ε = 0, 修正自动关闭 (区别于旧版恒定 δ = coeff*(|Ûn_rel|+ĉ))。
//   防止声速点/膨胀扇处特征值过零产生膨胀激波等非物理解。
// - 非物理保护 (SU2 roe.cpp): Roe 平均 ĉ² <= 0 时该面通量整体置零
//   返回 (启动瞬态保护; 正常计算经 MUSCL_WENO_Reconstruct.H 的
//   bad_recon 一阶回退后不应触发)。
// - 振动能量按被动标量迎风输运 (特征值 Ûn, 同样施加熵修正)。
//   ALE: 特征值用相对网格法向速度 (Ûn - u_mesh,n), 特征向量不变。
// 函数签名与 KurganovFlux / TadmorFlux / AUSMFlux 一致, 注册于同一
// fluxMap; face 值可由 KNP 机制 (reconstruct 格式) 或 MUSCL/WENO 重构
// 提供。
//------------------------------------------------------------------------------

#ifndef Roe_Flux_H
#define Roe_Flux_H

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// --- 熵修正配置 (由 readFluxScheme.H 从 fvSchemes 设置) ---
double roeEntropyFixCoeff = 1.0;

// --- Harten-Hyman 熵修正后的 |λ| ---
inline double roeFixAbsLambda(const double lambda, const double delta)
{
    const double absL = std::fabs(lambda);

    if (delta < 1.0e-30)
    {
        return absL;
    }

    return (absL < delta) ? (sqr(lambda) + sqr(delta))/(2.0*delta) : absL;
}


void RoeFlux
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
    const scalar magSf = mag(Sf);
    const vector n = Sf/magSf;

    // --- 面两侧总焓 (e 为比内能) ---
    const scalar H_pos = e_pos + p_pos/rho_pos + 0.5*magSqr(U_pos);
    const scalar H_neg = e_neg + p_neg/rho_neg + 0.5*magSqr(U_neg);

    // --- Roe 平均 (√ρ 加权) ---
    const scalar sL = Foam::sqrt(rho_pos);
    const scalar sR = Foam::sqrt(rho_neg);
    const scalar wL = sL/(sL + sR);
    const scalar wR = 1.0 - wL;

    const scalar rhoHat = sL*sR;
    const vector UHat = wL*U_pos + wR*U_neg;
    const scalar HHat = wL*H_pos + wR*H_neg;

    // 比振动能 (单温 0; 单振动模式 evk[0]; 多组分 Σ Yk*evk)
    scalar evHat_pos = 0.0;
    scalar evHat_neg = 0.0;
    if (!downgradeSingleT)
    {
        if (downgradeSingleTv)
        {
            evHat_pos = evk_pos[0];
            evHat_neg = evk_neg[0];
        }
        else
        {
            forAll(pD_pos, speciei)
            {
                evHat_pos += (pD_pos[speciei]/rho_pos)*evk_pos[speciei];
                evHat_neg += (pD_neg[speciei]/rho_neg)*evk_neg[speciei];
            }
        }
    }
    const scalar evHat = wL*evHat_pos + wR*evHat_neg;

    // 冻结比热比 (由 c² = γp/ρ 反解), √ρ 加权
    const scalar c_pos = cSf_pos/magSf;
    const scalar c_neg = cSf_neg/magSf;
    const scalar gamma_pos = rho_pos*sqr(c_pos)/max(p_pos, SMALL);
    const scalar gamma_neg = rho_neg*sqr(c_neg)/max(p_neg, SMALL);
    const scalar gammaHat = wL*gamma_pos + wR*gamma_neg;

    // --- 法向速度 (ALE: 相对网格) ---
    const scalar unHat = UHat & n;
    const scalar mesh_un = meshMoving ? mesh_phi/magSf : 0.0;
    const scalar unRel = unHat - mesh_un;
    const scalar unL = (U_pos & n) - mesh_un;
    const scalar unR = (U_neg & n) - mesh_un;

    // --- 声速 (两温近似, 见文件头) ---
    // 非物理保护 (SU2 roe.cpp): ĉ² <= 0 时不用 SMALL 截断 (会把波强
    // 放大至 ~1/SMALL 量级), 而是将该面通量整体置零返回。
    const scalar cHat2 =
        (gammaHat - 1.0)*(HHat - 0.5*magSqr(UHat) - evHat);

    if (cHat2 <= 0.0)
    {
        rhoFlux = 0.0;
        rhoUFlux = vector::zero;
        rhoEFlux = 0.0;

        if (downgradeSingleT)
        {
            rhoEvkFlux.clear();
        }
        else
        {
            forAll(rhoEvkFlux, speciei)
            {
                rhoEvkFlux[speciei] = 0.0;
            }
        }

        // CFL 用两侧最大波速 (Roe 平均不可用, 退回侧值)
        amaxSf = (max(mag(unL), mag(unR)) + max(c_pos, c_neg))*magSf;
        aByU = UHat;
        return;
    }

    const scalar cHat = Foam::sqrt(cHat2);

    // --- 跳跃量 ---
    const scalar drho = rho_neg - rho_pos;
    const vector dU = U_neg - U_pos;
    const scalar dun = dU & n;
    const scalar dp = p_neg - p_pos;
    const vector dUt = dU - dun*n;   // 切向速度跳变

    // --- 波强度 ---
    const scalar alpha1 = (dp - rhoHat*cHat*dun)/(2.0*cHat2);
    const scalar alpha2 = drho - dp/cHat2;
    const scalar alpha3 = (dp + rhoHat*cHat*dun)/(2.0*cHat2);

    // --- 熵修正后的特征值 |λ| (Harten-Hyman 逐波带宽, SU2 NEMO roe.cpp) ---
    // ε_k = EntropyCoeff*4*max(0, max(λ̂_k - λ_k,L, λ_k,R - λ̂_k)),
    // λ_k,L/R 为该波在面两侧状态的特征值 (网格相对速度)。
    const scalar epsU = roeEntropyFixCoeff*4.0*max
    (
        0.0,
        max(unRel - unL, unR - unRel)
    );
    const scalar epsM = roeEntropyFixCoeff*4.0*max
    (
        0.0,
        max((unRel - cHat) - (unL - c_pos), (unR - c_neg) - (unRel - cHat))
    );
    const scalar epsP = roeEntropyFixCoeff*4.0*max
    (
        0.0,
        max((unRel + cHat) - (unL + c_pos), (unR + c_neg) - (unRel + cHat))
    );

    const scalar l1 = roeFixAbsLambda(unRel - cHat, epsM);
    const scalar l2 = roeFixAbsLambda(unRel, epsU);
    const scalar l3 = roeFixAbsLambda(unRel + cHat, epsP);

    // --- 耗散项 (单位面积) ---
    const scalar Drho = l1*alpha1 + l2*alpha2 + l3*alpha3;

    const vector DU =
        l1*alpha1*(UHat - cHat*n)
      + l2*(alpha2*UHat + rhoHat*dUt)
      + l3*alpha3*(UHat + cHat*n);

    // 剪切波对能量通量的贡献: ρ̂(Û·Δu - Ûn·Δu_n) = ρ̂·Û_t·Δu_t
    const scalar shearE = rhoHat*((UHat & dU) - unHat*dun);
    const scalar DE =
        l1*alpha1*(HHat - unHat*cHat)
      + l2*(alpha2*0.5*magSqr(UHat) + shearE)
      + l3*alpha3*(HHat + unHat*cHat);

    // --- 对流通量 (相对网格速度) ---
    scalar phiv_pos = U_pos & Sf;
    scalar phiv_neg = U_neg & Sf;
    if (meshMoving)
    {
        phiv_pos -= mesh_phi;
        phiv_neg -= mesh_phi;
    }

    rhoFlux = 0.5*(phiv_pos*rho_pos + phiv_neg*rho_neg) - 0.5*magSf*Drho;

    rhoUFlux =
        0.5*(phiv_pos*rho_pos*U_pos + phiv_neg*rho_neg*U_neg)
      + 0.5*(p_pos + p_neg)*Sf
      - 0.5*magSf*DU;

    rhoEFlux =
        0.5*(phiv_pos*rho_pos*H_pos + phiv_neg*rho_neg*H_neg)
      - 0.5*magSf*DE;
    if (meshMoving)
    {
        rhoEFlux += mesh_phi*0.5*(p_pos + p_neg);
    }

    // --- 振动能量通量 (被动标量迎风, 特征值 unRel 施加熵修正) ---
    if (downgradeSingleT)
    {
        rhoEvkFlux.clear();
    }
    else if (downgradeSingleTv)
    {
        const scalar evcL = rho_pos*evk_pos[0];
        const scalar evcR = rho_neg*evk_neg[0];
        rhoEvkFlux[0] = 0.5*(phiv_pos*evcL + phiv_neg*evcR)
            - 0.5*l2*magSf*(evcR - evcL);
    }
    else
    {
        forAll(rhoEvkFlux, speciei)
        {
            const scalar evcL = pD_pos[speciei]*evk_pos[speciei];
            const scalar evcR = pD_neg[speciei]*evk_neg[speciei];
            rhoEvkFlux[speciei] = 0.5*(phiv_pos*evcL + phiv_neg*evcR)
                - 0.5*l2*magSf*(evcR - evcL);
        }
    }

    // --- CFL 用最大波速 / 粘性功用速度 ---
    amaxSf = (mag(unRel) + cHat)*magSf;
    aByU = UHat;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
