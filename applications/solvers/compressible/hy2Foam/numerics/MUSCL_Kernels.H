//------------------------------------------------------------------------------
// MUSCL-3 核函数 (κ-格式 + 总坡度限制器, SU2 对齐)
//
// 面 L/R 状态重构，用于 MUSCL_WENO_Reconstruct.H 中覆写 pos/neg 面场。
//
// 数学结构 (SU2 CEulerSolver::Upwind_Residual + CLimiterDetails.hpp):
//   κ 混合半坡度 (SU2 umusclProjection 的均匀网格等价形式):
//     P_L = 1/2*[(1-κ)*Δ_- + (1+κ)*Δ_+]
//     P_R = 1/2*[(1-κ)*Δ_2  + (1+κ)*Δ_+]
//   其中 Δ_- = u_C - u_{C-1}, Δ_+ = u_{C+1} - u_C (棱边差 δ),
//   Δ_2 = u_{C+2} - u_{C+1} (下游侧迎风差)。
//   重构 (SU2: V_L = V_i + 0.5*lim*Project_Grad_i):
//     u^L = u_C     + 0.5*L(P_L, Δ_+)
//     u^R = u_{C+1} - 0.5*L(P_R, Δ_+)
//   限制器作用于 [整个 κ 混合坡度] (旧版仅乘 (1+κ)Δ_+ 项, 极值肩部
//   残留 (1-κ)|Δ_-|/4 过冲, 且 κ < 1 时局部单调性可被破坏, 见 README)。
//
// κ 由 fvSchemes 字段 "MUSCLKappa" 指定 (默认 1/3):
//   κ =  1/3  3 阶 (标准 MUSCL 值)
//   κ =  0    Fromm 格式 (2 阶)
//   κ = -1    2 阶全迎风
//   κ =  1    2 阶中心
//
// 限制器由 fvSchemes 字段 "MUSCLLimiter" 指定 (默认 vanAlbada),
// 以坡度空间 (P, δ) 表示:
//   vanAlbada (SU2 vanAlbadaFunction, ε = 1e-16):
//     L = P*δ*(δ+P)/(P² + δ² + ε)
//   其余三个为同一 (P, δ) 槽位的经典坡度限制器:
//     minmod    L = minmod(P, δ)
//     vanLeer   L = (Pδ > 0) ? 2Pδ/(P+δ) : 0        (调和平均)
//     superbee  L = ±max(0, min(2|P|,|δ|), min(|P|,2|δ|)) (同号时)
// 性质 (四者一致):
//   - 光滑区 (P = δ): L = δ, 重构退化为无限制 κ-MUSCL
//     (κ = 1/3 时 3 阶), 与旧版光滑区结果完全一致;
//   - |0.5*L| <= |δ|: 重构值不越出棱边两端点值区间;
//   - 极值肩部 (δ → 0): L → 0, u^L → u_C, 无过冲。
//
// 提供标量和矢量版本。矢量版本逐分量施加 (SU2 对各原始变量亦
// 逐分量限制)。
//------------------------------------------------------------------------------

#ifndef MUSCL_Kernels_H
#define MUSCL_Kernels_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// --- 限制器与 κ 配置 (由 readFluxScheme.H 从 fvSchemes 设置) ---
// musclLimiterType: 0 = minmod, 1 = vanLeer, 2 = superbee, 3 = vanAlbada
int musclLimiterType = 3;
double musclKappa = 1.0/3.0;

namespace Foam
{

// --- 受限量坡度 L(P, δ) ---
// P: κ 混合半坡度; δ: 棱边差 (u_{C+1} - u_C), L/R 状态共用同一 δ
// (与 SU2 vanAlbadaFunction(Project_Grad, V_ij, EPS) 的参数约定一致)。
inline scalar musclLimitedSlope(const scalar P, const scalar delta)
{
    switch (musclLimiterType)
    {
        case 0: // minmod
        {
            return 0.5*(sign(P) + sign(delta))*min(mag(P), mag(delta));
        }
        case 1: // van Leer (调和平均)
        {
            return (P*delta > 0.0) ? 2.0*P*delta/(P + delta) : 0.0;
        }
        case 2: // superbee
        {
            if (P*delta <= 0.0)
            {
                return 0.0;
            }
            return sign(delta)
                *max(min(2.0*mag(P), mag(delta)), min(mag(P), 2.0*mag(delta)));
        }
        default: // van Albada (SU2 vanAlbadaFunction, EPS = 1e-16)
        {
            const scalar eps = 1.0e-16;
            return P*delta*(delta + P)/(sqr(P) + sqr(delta) + eps);
        }
    }
}


// --- MUSCL 左状态 (面 i+1/2, 上游侧) ---
// u^L = u_C + 0.5*L(P_L, Δ_+),  P_L = 1/2*[(1-κ)Δ_- + (1+κ)Δ_+]
inline scalar muscl3FaceL
(
    const scalar uM1, const scalar uC, const scalar uP1
)
{
    const scalar dM = uC - uM1;   // u_C - u_{C-1}
    const scalar dP = uP1 - uC;   // u_{C+1} - u_C  (棱边差 δ)
    const scalar P = 0.5*((1.0 - musclKappa)*dM + (1.0 + musclKappa)*dP);
    return uC + 0.5*musclLimitedSlope(P, dP);
}

inline vector muscl3FaceL
(
    const vector& uM1, const vector& uC, const vector& uP1
)
{
    vector uL(uC);
    for (label cmp = 0; cmp < 3; ++cmp)
    {
        uL[cmp] = muscl3FaceL(uM1[cmp], uC[cmp], uP1[cmp]);
    }
    return uL;
}


// --- MUSCL 右状态 (面 i+1/2, 下游侧) ---
// u^R = u_{C+1} - 0.5*L(P_R, Δ_+),  P_R = 1/2*[(1-κ)Δ_2 + (1+κ)Δ_+]
// Δ_+ = u_{C+1} - u_C 为棱边差, 与 L 状态共用 (SU2 亦如此);
// muscl3FaceR(a,b,c) 与 muscl3FaceL(c,b,a) 网格反射镜像。
inline scalar muscl3FaceR
(
    const scalar uC, const scalar uP1, const scalar uP2
)
{
    const scalar d1 = uP1 - uC;   // 棱边差 δ (界面侧)
    const scalar d2 = uP2 - uP1;  // u_{C+2} - u_{C+1} (迎风侧)
    const scalar P = 0.5*((1.0 - musclKappa)*d2 + (1.0 + musclKappa)*d1);
    return uP1 - 0.5*musclLimitedSlope(P, d1);
}

inline vector muscl3FaceR
(
    const vector& uC, const vector& uP1, const vector& uP2
)
{
    vector uR(uP1);
    for (label cmp = 0; cmp < 3; ++cmp)
    {
        uR[cmp] = muscl3FaceR(uC[cmp], uP1[cmp], uP2[cmp]);
    }
    return uR;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
