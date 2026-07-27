//------------------------------------------------------------------------------
// WENO-Z 核函数 (Borges et al., JCP 2008)
//
// WENO-Z-5: 5 单元模板 → 3 个 3 单元候选模板, γ=(0.1, 0.6, 0.3)
// WENO-Z-7: 7 单元模板 → 4 个 4 单元候选模板, γ=(1/35, 12/35, 18/35, 4/35)
//
// 标量和矢量版本。矢量版本逐分量应用 WENO 权重。
//------------------------------------------------------------------------------

#ifndef WENOZ_Kernels_H
#define WENOZ_Kernels_H

#include "scalar.H"
#include "vector.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// ====================================================================
// WENO-Z-5
// ====================================================================

// --- 左状态: u^L_{i+1/2} ---
// 5 点模板: {uM2, uM1, uC, uP1, uP2}
// 3 个候选模板的插值
inline scalar wenoZ5_faceL
(
    const scalar uM2, const scalar uM1, const scalar uC,
    const scalar uP1, const scalar uP2
)
{
    const scalar eps = 1e-40;

    // 候选模板插值 p_k
    const scalar p0 = ( 2.0*uM2 -  7.0*uM1 + 11.0*uC )/6.0;
    const scalar p1 = (      -1.0*uM1 +  5.0*uC +  2.0*uP1)/6.0;
    const scalar p2 = (                   2.0*uC +  5.0*uP1 - uP2)/6.0;

    // 光滑度指示器 β_k
    const scalar b0 = (13.0/12.0)*sqr(uM2 - 2.0*uM1 + uC)
                    + (1.0/4.0)*sqr(uM2 - 4.0*uM1 + 3.0*uC);
    const scalar b1 = (13.0/12.0)*sqr(uM1 - 2.0*uC + uP1)
                    + (1.0/4.0)*sqr(uM1 - uP1);
    const scalar b2 = (13.0/12.0)*sqr(uC - 2.0*uP1 + uP2)
                    + (1.0/4.0)*sqr(3.0*uC - 4.0*uP1 + uP2);

    // 全局光滑度参考 τ₅ = |β₀ - β₂|
    const scalar tau5 = mag(b0 - b2);

    // 线性权重
    const scalar g0 = 0.1;
    const scalar g1 = 0.6;
    const scalar g2 = 0.3;

    // WENO-Z 权重: α_k = γ_k * (1 + τ₅/(β_k + ε))
    const scalar a0 = g0*(1.0 + tau5/(b0 + eps));
    const scalar a1 = g1*(1.0 + tau5/(b1 + eps));
    const scalar a2 = g2*(1.0 + tau5/(b2 + eps));
    const scalar aSum = a0 + a1 + a2;

    // 归一化权重
    const scalar w0 = a0/aSum;
    const scalar w1 = a1/aSum;
    const scalar w2 = a2/aSum;

    return w0*p0 + w1*p1 + w2*p2;
}


// --- 右状态: u^R_{i+1/2} ---
// 5 点模板: {uM1, uC, uP1, uP2, uP3} (注意右偏)
inline scalar wenoZ5_faceR
(
    const scalar uM1, const scalar uC, const scalar uP1,
    const scalar uP2, const scalar uP3
)
{
    const scalar eps = 1e-40;

    // 候选模板插值 (右偏)
    // S0^R = {uP1, uP2, uP3}:
    const scalar p0 = ( 2.0*uP3 -  7.0*uP2 + 11.0*uP1)/6.0;
    // S1^R = {uC, uP1, uP2}:
    const scalar p1 = (      -1.0*uP2 +  5.0*uP1 +  2.0*uC)/6.0;
    // S2^R = {uM1, uC, uP1}:
    const scalar p2 = (                   2.0*uP1 +  5.0*uC - uM1)/6.0;

    // 光滑度指示器 (右偏)
    const scalar b0 = (13.0/12.0)*sqr(uP3 - 2.0*uP2 + uP1)
                    + (1.0/4.0)*sqr(uP3 - 4.0*uP2 + 3.0*uP1);
    const scalar b1 = (13.0/12.0)*sqr(uP2 - 2.0*uP1 + uC)
                    + (1.0/4.0)*sqr(uP2 - uC);
    const scalar b2 = (13.0/12.0)*sqr(uP1 - 2.0*uC + uM1)
                    + (1.0/4.0)*sqr(3.0*uP1 - 4.0*uC + uM1);

    const scalar tau5 = mag(b0 - b2);

    const scalar g0 = 0.1;
    const scalar g1 = 0.6;
    const scalar g2 = 0.3;

    const scalar a0 = g0*(1.0 + tau5/(b0 + eps));
    const scalar a1 = g1*(1.0 + tau5/(b1 + eps));
    const scalar a2 = g2*(1.0 + tau5/(b2 + eps));
    const scalar aSum = a0 + a1 + a2;

    return (a0*p0 + a1*p1 + a2*p2)/aSum;
}


// ====================================================================
// WENO-Z-7
// ====================================================================

// --- 4 阶子模板的插值系数 (左状态, 面 i+1/2) ---
// v0: S0 = {i-3, i-2, i-1, i}  → u^L_{i+1/2} = (-3u_{i-3} + 13u_{i-2} - 23u_{i-1} + 25u_i)/12
// v1: S1 = {i-2, i-1, i, i+1}  → u^L_{i+1/2} = (u_{i-2} - 5u_{i-1} + 13u_i + 3u_{i+1})/12
// v2: S2 = {i-1, i, i+1, i+2}  → u^L_{i+1/2} = (-u_{i-1} + 7u_i + 7u_{i+1} - u_{i+2})/12
// v3: S3 = {i, i+1, i+2, i+3}  → u^L_{i+1/2} = (3u_i + 13u_{i+1} - 5u_{i+2} + u_{i+3})/12

// --- 4-cell smoothness indicator (Jiang-Shu 1996) ---
// β = c1*(v1²+v2²+v3²) + c2*(d1²+d2²) + c3*t²
//   where v_i are 1st undivided diffs, d_i 2nd diffs, t the 3rd diff
inline scalar beta4
(
    const scalar f0, const scalar f1, const scalar f2, const scalar f3
)
{
    // 1st undivided differences
    const scalar v1 = f1 - f0;
    const scalar v2 = f2 - f1;
    const scalar v3 = f3 - f2;

    // 2nd undivided differences
    const scalar d1 = f0 - 2.0*f1 + f2;
    const scalar d2 = f1 - 2.0*f2 + f3;

    // 3rd undivided difference
    const scalar t = -f0 + 3.0*f1 - 3.0*f2 + f3;

    // 系数通过 Taylor 展开确定 (对应 4 阶多项式的光滑度度量)
    // 这些系数确保 β 在光滑区域 ∝ O(Δx²)
    return (v1*(547.0*v1 - 2522.0*v2 + 3441.0*v3)
          + v2*(3043.0*v2 - 8534.0*v3)
          + v3*(6243.0*v3)
          + d1*(547.0*d1 - 2522.0*d2)
          + d2*(3043.0*d2)
          + 547.0*sqr(t)) / 240.0;
}


// --- 左状态: u^L_{i+1/2}, 7 点模板 ---
// 输入: {uM3, uM2, uM1, uC, uP1, uP2, uP3}
inline scalar wenoZ7_faceL
(
    const scalar uM3, const scalar uM2, const scalar uM1, const scalar uC,
    const scalar uP1, const scalar uP2, const scalar uP3
)
{
    const scalar eps = 1e-40;

    // 候选模板插值
    const scalar v0 = (-3.0*uM3 + 13.0*uM2 - 23.0*uM1 + 25.0*uC)/12.0;
    const scalar v1 = (       1.0*uM2 -  5.0*uM1 + 13.0*uC +  3.0*uP1)/12.0;
    const scalar v2 = (                 -1.0*uM1 +  7.0*uC +  7.0*uP1 - uP2)/12.0;
    const scalar v3 = (                             3.0*uC + 13.0*uP1 - 5.0*uP2 + uP3)/12.0;

    // 光滑度指示器
    const scalar b0 = beta4(uM3, uM2, uM1, uC);
    const scalar b1 = beta4(uM2, uM1, uC, uP1);
    const scalar b2 = beta4(uM1, uC, uP1, uP2);
    const scalar b3 = beta4(uC, uP1, uP2, uP3);

    // 全局光滑度 τ₇ = |β₀ - 3β₁ + 3β₂ - β₃|
    const scalar tau7 = mag(b0 - 3.0*b1 + 3.0*b2 - b3);

    // 线性权重
    const scalar g0 =  1.0/35.0;
    const scalar g1 = 12.0/35.0;
    const scalar g2 = 18.0/35.0;
    const scalar g3 =  4.0/35.0;

    // WENO-Z 权重
    const scalar a0 = g0*(1.0 + tau7/(b0 + eps));
    const scalar a1 = g1*(1.0 + tau7/(b1 + eps));
    const scalar a2 = g2*(1.0 + tau7/(b2 + eps));
    const scalar a3 = g3*(1.0 + tau7/(b3 + eps));
    const scalar aSum = a0 + a1 + a2 + a3;

    return (a0*v0 + a1*v1 + a2*v2 + a3*v3)/aSum;
}


// --- 右状态: u^R_{i+1/2}, 7 点模板 ---
// 输入: {uM2, uM1, uC, uP1, uP2, uP3, uP4}
inline scalar wenoZ7_faceR
(
    const scalar uM2, const scalar uM1, const scalar uC, const scalar uP1,
    const scalar uP2, const scalar uP3, const scalar uP4
)
{
    const scalar eps = 1e-40;

    // 候选模板插值 (右偏)
    // S0^R = {i+4, i+3, i+2, i+1}:
    const scalar v0 = (-3.0*uP4 + 13.0*uP3 - 23.0*uP2 + 25.0*uP1)/12.0;
    // S1^R = {i+3, i+2, i+1, i}:
    const scalar v1 = (       1.0*uP3 -  5.0*uP2 + 13.0*uP1 +  3.0*uC)/12.0;
    // S2^R = {i+2, i+1, i, i-1}:
    const scalar v2 = (                 -1.0*uP2 +  7.0*uP1 +  7.0*uC - uM1)/12.0;
    // S3^R = {i+1, i, i-1, i-2}:
    const scalar v3 = (                             3.0*uP1 + 13.0*uC - 5.0*uM1 + uM2)/12.0;

    // 光滑度指示器
    const scalar b0 = beta4(uP4, uP3, uP2, uP1);
    const scalar b1 = beta4(uP3, uP2, uP1, uC);
    const scalar b2 = beta4(uP2, uP1, uC, uM1);
    const scalar b3 = beta4(uP1, uC, uM1, uM2);

    const scalar tau7 = mag(b0 - 3.0*b1 + 3.0*b2 - b3);

    const scalar g0 =  1.0/35.0;
    const scalar g1 = 12.0/35.0;
    const scalar g2 = 18.0/35.0;
    const scalar g3 =  4.0/35.0;

    const scalar a0 = g0*(1.0 + tau7/(b0 + eps));
    const scalar a1 = g1*(1.0 + tau7/(b1 + eps));
    const scalar a2 = g2*(1.0 + tau7/(b2 + eps));
    const scalar a3 = g3*(1.0 + tau7/(b3 + eps));
    const scalar aSum = a0 + a1 + a2 + a3;

    return (a0*v0 + a1*v1 + a2*v2 + a3*v3)/aSum;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
