// --- 处理器边界 ghost 层构建 (并行 MUSCL/WENO 模板) ---
// 初始化时执行一次, 紧随 structuredStencilBuild.H 之后在同一个
// if (fluxScheme == "AUSM"/"Roe") 块内被 readFluxScheme.H include。
//
// 目标: 让 cell_±_x/y/z 方向链在处理器边界处延伸到远端进程的单元
// (ghost 单元, 最多 ghostDepth 层), 使高阶模板在分区边界不再退化为
// 常值外推。ghost 单元的场值由 ghostFieldExchange.H 每时间步交换。
//
// ghost 单元编号 (全局): nLocalCells + ghostStart[ppI] + face*ghostDepth
// + (d-1), 其中 d = 1..ghostDepth 为距处理器边界的层数。
// processor patch 的面序与耦合对端一一对应 (OpenFOAM 保证), 因此
// 发送缓冲按 (face, d) 排序即可与接收方的 ghost 编号对齐。
//
// 变量 (由 readFluxScheme.H 声明):
//   ghostDepth, nGhostCells, procPatchIDs, ghostSendCells, ghostStart

// ghost 深度: 模板最远到达 N+3 (WENOZ7) -> 3 层; WENOZ5 -> 2; MUSCL -> 1
ghostDepth = (convSchemes == "WENOZ7") ? 3
           : (convSchemes == "WENOZ5") ? 2 : 1;

const polyBoundaryMesh& patchesG = mesh.boundaryMesh();

// 收集参与交换的 processor patch (串行或无邻接 -> 空表, nGhostCells = 0)
label nProcPatchesG = 0;
forAll(patchesG, patchi)
{
    if (isA<processorPolyPatch>(patchesG[patchi]) and patchesG[patchi].size() > 0)
    {
        nProcPatchesG++;
    }
}
procPatchIDs.setSize(nProcPatchesG);
ghostSendCells.setSize(nProcPatchesG);
ghostStart.setSize(nProcPatchesG);
{
    label ppI = 0;
    forAll(patchesG, patchi)
    {
        if
        (
            isA<processorPolyPatch>(patchesG[patchi])
         and patchesG[patchi].size() > 0
        )
        {
            procPatchIDs[ppI++] = patchi;
        }
    }
}

const surfaceVectorField& SfG = mesh.Sf();
const surfaceScalarField& magSfG = mesh.magSf();

// 每面外向轴标记 (供第二遍链延伸复用): 0=非轴对齐, 其余同 axisFaces 约定
labelListList outwardTagsG(nProcPatchesG);

// 非轴对齐处理器面计数 (不延伸链, 保留常值外推行为)
label nNonAxisProcFaces = 0;

// --- 第一遍: 计算发送走链 (此时尚未扩展链, 走向只经本地单元) ---
// 发送给对端的单元 = 本侧从边界面 owner 沿域内方向走 0..ghostDepth-1 步
// (即对端进程的 ghost 层)。分区内向走链遇 -1 时重复最后有效单元,
// 与串行物理边界的常值外推语义一致。
nGhostCells = 0;

forAll(procPatchIDs, ppI)
{
    const label patchi = procPatchIDs[ppI];
    const labelUList& faceCellsG = patchesG[patchi].faceCells();
    const vectorField& SfbG = SfG.boundaryField()[patchi];
    const scalarField& magSfbG = magSfG.boundaryField()[patchi];
    const label nFacesG = faceCellsG.size();

    outwardTagsG[ppI].setSize(nFacesG, 0);
    ghostStart[ppI] = nGhostCells;
    nGhostCells += nFacesG*ghostDepth;
    ghostSendCells[ppI].setSize(nFacesG*ghostDepth);

    for (label facei = 0; facei < nFacesG; ++facei)
    {
        const label c = faceCellsG[facei];

        // 轴向判定 (与 structuredStencilBuild 相同阈值 0.9)
        // Sf 指向 owner -> 远端, 故法向符号即外向方向
        label outwardTag = 0;
        if (magSfbG[facei] > VSMALL)
        {
            const vector nG = SfbG[facei]/magSfbG[facei];
            const scalar dotXG = mag(nG & vector(1,0,0));
            const scalar dotYG = mag(nG & vector(0,1,0));
            const scalar dotZG = mag(nG & vector(0,0,1));

            if (dotXG > 0.9)
            {
                outwardTag = ((nG & vector(1,0,0)) > 0) ? 1 : -1;
            }
            else if (dotYG > 0.9)
            {
                outwardTag = ((nG & vector(0,1,0)) > 0) ? 2 : -2;
            }
            else if (dotZG > 0.9)
            {
                outwardTag = ((nG & vector(0,0,1)) > 0) ? 3 : -3;
            }
        }
        outwardTagsG[ppI][facei] = outwardTag;

        if (outwardTag == 0)
        {
            // 非轴对齐处理器面: 值不会被引用 (链不延伸), 填充 owner 保序
            nNonAxisProcFaces++;
            for (label d = 1; d <= ghostDepth; ++d)
            {
                ghostSendCells[ppI][facei*ghostDepth + d - 1] = c;
            }
            continue;
        }

        // 发送走链: getDirLists(outwardTag) 的 plusDir 沿外向,
        // offset 取负即向域内走; d=1 (offset 0) 为 owner 本身
        const labelList* plusDirG = nullptr;
        const labelList* minusDirG = nullptr;
        getDirLists
        (
            outwardTag,
            cell_plus_x, cell_minus_x,
            cell_plus_y, cell_minus_y,
            cell_plus_z, cell_minus_z,
            plusDirG, minusDirG
        );

        for (label d = 1; d <= ghostDepth; ++d)
        {
            ghostSendCells[ppI][facei*ghostDepth + d - 1] =
                ghostStencilCell(c, -(d-1), *plusDirG, *minusDirG);
        }
    }
}

// --- 扩展六条方向链: 追加 ghost 区 (先全部置 -1) ---
const label nLocalCellsG = mesh.nCells();
const label nTotalCellsG = nLocalCellsG + nGhostCells;

cell_plus_x.setSize(nTotalCellsG);
cell_minus_x.setSize(nTotalCellsG);
cell_plus_y.setSize(nTotalCellsG);
cell_minus_y.setSize(nTotalCellsG);
cell_plus_z.setSize(nTotalCellsG);
cell_minus_z.setSize(nTotalCellsG);

for (label c = nLocalCellsG; c < nTotalCellsG; ++c)
{
    cell_plus_x[c]  = -1;
    cell_minus_x[c] = -1;
    cell_plus_y[c]  = -1;
    cell_minus_y[c] = -1;
    cell_plus_z[c]  = -1;
    cell_minus_z[c] = -1;
}

// --- 第二遍: 链延伸 --- 边界单元外向槽指向 ghost 链
// ghost g_d --外向--> g_{d+1} (d < ghostDepth), 最外层保持 -1;
// 同时补内向: g_1 的内侧邻为边界单元, g_{d+1} 的内侧邻为 g_d
// (使走链在 ghost 区内自洽)。
label nChainClashes = 0;

forAll(procPatchIDs, ppI)
{
    const label patchi = procPatchIDs[ppI];
    const labelUList& faceCellsG = patchesG[patchi].faceCells();
    const labelList& tagsG = outwardTagsG[ppI];

    forAll(faceCellsG, facei)
    {
        const label outwardTag = tagsG[facei];
        if (outwardTag == 0) continue;

        const label c = faceCellsG[facei];
        const label gBase = nLocalCellsG + ghostStart[ppI] + facei*ghostDepth;

        // 外向/内向列表 (按外向符号选 cell_±_k)
        labelList* outListG = nullptr;
        labelList* inListG = nullptr;
        switch (outwardTag)
        {
            case  1: outListG = &cell_plus_x;  inListG = &cell_minus_x; break;
            case -1: outListG = &cell_minus_x; inListG = &cell_plus_x;  break;
            case  2: outListG = &cell_plus_y;  inListG = &cell_minus_y; break;
            case -2: outListG = &cell_minus_y; inListG = &cell_plus_y;  break;
            case  3: outListG = &cell_plus_z;  inListG = &cell_minus_z; break;
            case -3: outListG = &cell_minus_z; inListG = &cell_plus_z;  break;
            default: continue;
        }

        if ((*outListG)[c] != -1)
        {
            // 外向槽已被占用 (网格拓扑异常): 跳过, 保留原链
            nChainClashes++;
            continue;
        }

        (*outListG)[c] = gBase;
        (*inListG)[gBase] = c;
        for (label d = 1; d < ghostDepth; ++d)
        {
            (*outListG)[gBase + d - 1] = gBase + d;
            (*inListG)[gBase + d] = gBase + d - 1;
        }
    }
}

// --- 输出汇总 ---
reduce(nGhostCells, sumOp<label>());
Info<< "Processor ghost layers: " << nProcPatchesG << " patches"
    << ", depth " << ghostDepth
    << ", " << nGhostCells << " ghost cells"
    << (nGhostCells > 0 ? "" : " (serial or single-domain)")
    << endl;

if (nNonAxisProcFaces > 0)
{
    WarningIn
    (
        "ghostStencilBuild"
    )   << nNonAxisProcFaces << " processor-boundary faces are not "
        << "axis-aligned: their stencils keep constant extrapolation "
        << "(same behaviour as physical boundaries)."
        << endl;
}

if (nChainClashes > 0)
{
    WarningIn
    (
        "ghostStencilBuild"
    )   << nChainClashes << " processor-boundary cells had a non-empty "
        << "outward neighbour slot: ghost chain skipped for those faces "
        << "(stencil keeps constant extrapolation)."
        << endl;
}

// ************************************************************************* //
