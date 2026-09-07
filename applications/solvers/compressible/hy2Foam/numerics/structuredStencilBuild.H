// 初始化时执行一次，通过面法向检测轴对齐面 (blockMesh 六面体网格)
// 构建每个单元沿 ±x/±y/±z 方向的邻居列表。
//
// 变量 (由 readFluxScheme.H 声明):
//   axisFaces, nAxisFaces, isStructuredGrid
//   cell_plus_x, cell_minus_x, cell_plus_y, cell_minus_y, cell_plus_z, cell_minus_z

axisFaces.setSize(mesh.nInternalFaces(), 0);
nAxisFaces = 0;

cell_plus_x.setSize (mesh.nCells(), -1);
cell_minus_x.setSize(mesh.nCells(), -1);
cell_plus_y.setSize (mesh.nCells(), -1);
cell_minus_y.setSize(mesh.nCells(), -1);
cell_plus_z.setSize (mesh.nCells(), -1);
cell_minus_z.setSize(mesh.nCells(), -1);

const surfaceVectorField& Sf = mesh.Sf();
const vectorField& C = mesh.C();

forAll(mesh.owner(), facei)
{
    const label own  = mesh.owner()[facei];
    const label nei  = mesh.neighbour()[facei];
    const vector n   = Sf[facei]/mesh.magSf()[facei];
    const vector dir = (C[nei] - C[own]);
    const scalar magDir = mag(dir);

    if (magDir < VSMALL)
    {
        continue;
    }

    // 面法向与 x/y/z 轴的点积
    const scalar dotX = mag(n & vector(1,0,0));
    const scalar dotY = mag(n & vector(0,1,0));
    const scalar dotZ = mag(n & vector(0,0,1));

    // 阈值 0.9: blockMesh 六面体网格的内部面精确轴对齐
    if (dotX > 0.9)
    {
        const scalar proj = dir & vector(1,0,0);
        if (proj > 0)
        {
            cell_plus_x[own]  = nei;
            cell_minus_x[nei] = own;
            axisFaces[facei]  = 1;
        }
        else
        {
            cell_plus_x[nei]  = own;
            cell_minus_x[own] = nei;
            axisFaces[facei]  = -1;
        }
        nAxisFaces++;
    }
    else if (dotY > 0.9)
    {
        const scalar proj = dir & vector(0,1,0);
        if (proj > 0)
        {
            cell_plus_y[own]  = nei;
            cell_minus_y[nei] = own;
            axisFaces[facei]  = 2;
        }
        else
        {
            cell_plus_y[nei]  = own;
            cell_minus_y[own] = nei;
            axisFaces[facei]  = -2;
        }
        nAxisFaces++;
    }
    else if (dotZ > 0.9)
    {
        const scalar proj = dir & vector(0,0,1);
        if (proj > 0)
        {
            cell_plus_z[own]  = nei;
            cell_minus_z[nei] = own;
            axisFaces[facei]  = 3;
        }
        else
        {
            cell_plus_z[nei]  = own;
            cell_minus_z[own] = nei;
            axisFaces[facei]  = -3;
        }
        nAxisFaces++;
    }
}

// 输出结构化网格检测结果
reduce(nAxisFaces, sumOp<label>());
const label nInternalFaces = returnReduce(mesh.nInternalFaces(), sumOp<label>());

Info<< "Structured stencil: " << nAxisFaces << " / "
    << nInternalFaces << " axis-aligned internal faces ("
    << (nInternalFaces > 0 ? 100.0*nAxisFaces/nInternalFaces : 0.0)
    << " %)" << endl;

// 标记是否使用结构化模板 (要求全部内部面轴对齐)
isStructuredGrid = (nInternalFaces > 0)
                && (nAxisFaces == nInternalFaces);

// MUSCL/WENO 依赖轴对齐方向链模板: 非轴对齐面会静默回落到二阶 TVD 插值,
// 丢失高阶精度, 因此检测不通过时直接退出 (不再静默降级)
if (!isStructuredGrid)
{
    FatalErrorIn
    (
        "structuredStencilBuild"
    )   << "fluxScheme " << fluxScheme
        << " (convSchemes " << convSchemes << ") requires a fully "
        << "axis-aligned structured hexahedral mesh, but only "
        << nAxisFaces << " of " << nInternalFaces << " internal faces ("
        << (nInternalFaces > 0 ? 100.0*nAxisFaces/nInternalFaces : 0.0)
        << " %) are axis-aligned." << nl
        << "Non-axis-aligned faces would silently fall back to 2nd-order "
        << "TVD interpolation and lose the requested high-order accuracy." << nl
        << "Use fluxScheme Kurganov (works on any mesh) or regenerate the "
        << "mesh with axis-aligned blocks (e.g. blockMesh)."
        << abort(FatalError);
}
