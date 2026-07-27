#ifndef highOrderKernels_H
#define highOrderKernels_H

#include "scalar.H"
#include "vector.H"
#include "labelList.H"

#include "MUSCL_Kernels.H"
#include "WENOZ_Kernels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// --- 虚网格辅助: 沿方向链查找 offset 步后的单元 ---
// direction > 0: 沿 cell_plus_d 方向, direction < 0: 沿 cell_minus_d 方向
// 若遇到边界 (-1), 返回最后有效单元 (常量外推)
inline Foam::label ghostStencilCell
(
    const Foam::label startCelli,
    const Foam::label offset,
    const Foam::labelList& plusDir,
    const Foam::labelList& minusDir
)
{
    Foam::label cur = startCelli;

    if (offset > 0)
    {
        for (Foam::label s = 0; s < offset; ++s)
        {
            const Foam::label nxt = plusDir[cur];
            if (nxt == -1) break;
            cur = nxt;
        }
    }
    else if (offset < 0)
    {
        for (Foam::label s = 0; s < -offset; ++s)
        {
            const Foam::label nxt = minusDir[cur];
            if (nxt == -1) break;
            cur = nxt;
        }
    }
    return cur;
}


// --- 标量场虚网格模板取值 ---
inline Foam::scalar stencilScalar
(
    const Foam::scalarField& fld,
    const Foam::label startCelli,
    const Foam::label offset,
    const Foam::labelList& plusDir,
    const Foam::labelList& minusDir
)
{
    return fld[ghostStencilCell(startCelli, offset, plusDir, minusDir)];
}


// --- 矢量场虚网格模板取值 ---
inline Foam::vector stencilVector
(
    const Foam::vectorField& fld,
    const Foam::label startCelli,
    const Foam::label offset,
    const Foam::labelList& plusDir,
    const Foam::labelList& minusDir
)
{
    return fld[ghostStencilCell(startCelli, offset, plusDir, minusDir)];
}


// --- 根据 face 方向选取 ± 邻居列表 ---
inline void getDirLists
(
    const Foam::label faceDir,
    const Foam::labelList& cell_plus_x,
    const Foam::labelList& cell_minus_x,
    const Foam::labelList& cell_plus_y,
    const Foam::labelList& cell_minus_y,
    const Foam::labelList& cell_plus_z,
    const Foam::labelList& cell_minus_z,
    const Foam::labelList*& plusDir,
    const Foam::labelList*& minusDir
)
{
    switch (faceDir)
    {
        case  1: plusDir = &cell_plus_x;  minusDir = &cell_minus_x; break;
        case -1: plusDir = &cell_minus_x; minusDir = &cell_plus_x;  break;
        case  2: plusDir = &cell_plus_y;  minusDir = &cell_minus_y; break;
        case -2: plusDir = &cell_minus_y; minusDir = &cell_plus_y;  break;
        case  3: plusDir = &cell_plus_z;  minusDir = &cell_minus_z; break;
        case -3: plusDir = &cell_minus_z; minusDir = &cell_plus_z;  break;
        default: plusDir = nullptr; minusDir = nullptr; break;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
