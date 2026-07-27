//------------------------------------------------------------------------------
// MUSCL-3 核函数 (κ=1/3 + van Albada 限制器)
//
// 面 L/R 状态重构，用于 MUSCL_WENO_Reconstruct.H 中覆写 pos/neg 面场。
// 提供标量和矢量版本。
//------------------------------------------------------------------------------

#ifndef MUSCL_Kernels_H
#define MUSCL_Kernels_H

#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// --- van Albada 限制器 ---
// 标量版: r = dM/(dP+ε)
inline scalar vanAlbadaLimiter(const scalar dM, const scalar dP)
{
    const scalar eps = 1e-20;
    const scalar r = dM/(dP + eps);
    // φ(r) = (r² + r) / (r² + 1) = r*(r+1)/(r²+1)
    return r*(r + 1.0)/(sqr(r) + 1.0);
}

// 矢量版: r = (dM·dP) / (|dP|² + ε)
inline scalar vanAlbadaLimiter(const vector& dM, const vector& dP)
{
    const scalar eps = 1e-20;
    const scalar r = (dM & dP)/(magSqr(dP) + eps);
    return r*(r + 1.0)/(sqr(r) + 1.0);
}


// --- MUSCL-3 左状态 (面 i+1/2, 上游侧) ---
// uL = u_C + (1/6)*dM + (1/3)*φ*dP
inline scalar muscl3FaceL
(
    const scalar uM1, const scalar uC, const scalar uP1
)
{
    const scalar dM = uC - uM1;   // u_C - u_{C-1}
    const scalar dP = uP1 - uC;   // u_{C+1} - u_C
    const scalar phi = vanAlbadaLimiter(dM, dP);
    return uC + (1.0/6.0)*dM + (1.0/3.0)*phi*dP;
}

inline vector muscl3FaceL
(
    const vector& uM1, const vector& uC, const vector& uP1
)
{
    const vector dM = uC - uM1;
    const vector dP = uP1 - uC;
    const scalar phi = vanAlbadaLimiter(dM, dP);
    return uC + (1.0/6.0)*dM + (1.0/3.0)*phi*dP;
}


// --- MUSCL-3 右状态 (面 i+1/2, 下游侧) ---
// uR = u_{C+1} - (1/6)*dP - (1/3)*φ*dM
// 其中 dM,uP1 = u_{C+1} - u_C, dP,uR = u_{C+2} - u_{C+1}
// 即 uR = u_P1 - (1/6)*(u_P2-u_P1) - (1/3)*φ*(u_P1-u_C)
inline scalar muscl3FaceR
(
    const scalar uC, const scalar uP1, const scalar uP2
)
{
    const scalar dM = uP1 - uC;   // u_{C+1} - u_C (C is N from R perspective)
    const scalar dP = uP2 - uP1;  // u_{C+2} - u_{C+1}
    const scalar phi = vanAlbadaLimiter(dM, dP);
    return uP1 - (1.0/6.0)*dP - (1.0/3.0)*phi*dM;
}

inline vector muscl3FaceR
(
    const vector& uC, const vector& uP1, const vector& uP2
)
{
    const vector dM = uP1 - uC;
    const vector dP = uP2 - uP1;
    const scalar phi = vanAlbadaLimiter(dM, dP);
    return uP1 - (1.0/6.0)*dP - (1.0/3.0)*phi*dM;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
