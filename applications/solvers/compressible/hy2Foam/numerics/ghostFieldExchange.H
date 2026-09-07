//------------------------------------------------------------------------------
// 处理器边界 ghost 层交换 (每时间步, 在 MUSCL/WENO 面重构之前) — 执行代码
//
// 由 MUSCL_WENO_Reconstruct.H 在其 fluxScheme = AUSM/Roe 块顶部 include。
// 组装 "本地 + ghost" 扩展单元场 (rho, rhoU, rPsi, e, evk, pD), 并与
// 相邻进程交换 ghost 层单元值, 使高阶模板在处理器边界处取得真实的
// 远端单元值 (而非常值外推)。
//
// 串行 (nGhostCells == 0) 时扩展场即原始单元场的拷贝, 走链与数值路径
// 完全不变 (与未加本文件时 bit 级一致)。
//
// ghost 映射与扩展后的方向链由 ghostStencilBuild.H 在初始化时构建;
// 变量 nGhostCells, procPatchIDs, ghostSendCells, ghostStart 于
// readFluxScheme.H 声明。交换采用 PstreamBuffers 非阻塞模式;
// 发送缓冲布局: [face*ghostDepth + (d-1)]*nVals + 槽位, 值槽依次为
// rho, rhoU(3 分量), rPsi, e, evk[nEv], pD[nPd] — 收发双方布局一致
// (processor patch 面序一一对应; 降解模式与物种数全场一致)。
//------------------------------------------------------------------------------

const label nLocalCellsX = mesh.nCells();
const label nTotalCellsX = nLocalCellsX + nGhostCells;

// --- 扩展场: 本地部分拷贝 primitiveField, ghost 部分由交换写入 ---
// 定义为拥有所有权的局部场 (而非引用), 以便串行/并行统一使用;
// nGhostCells == 0 时 setSize 为空操作, 即原场拷贝
scalarField rhoI(rho.primitiveField());
vectorField rhoUI(rhoU.primitiveField());
scalarField rPsiI(rPsi.primitiveField());
scalarField eI(e.primitiveField());

rhoI.setSize(nTotalCellsX);
rhoUI.setSize(nTotalCellsX);
rPsiI.setSize(nTotalCellsX);
eI.setSize(nTotalCellsX);

// 振动能扩展场: 单振动模式为 ev, 多物种模式为 hevk (与
// MUSCL_WENO_Reconstruct.H 的 evk_pos 索引约定一致)
PtrList<scalarField> evkExt(evk_pos.size());
PtrList<scalarField> pDExt(pD_pos.size());

forAll(evkExt, speciei)
{
    if (downgradeSingleTv)
    {
        evkExt.set(speciei, new scalarField(ev.primitiveField()));
    }
    else
    {
        evkExt.set
        (
            speciei,
            new scalarField(hevk[speciei].primitiveField())
        );
    }
    evkExt[speciei].setSize(nTotalCellsX);
}

forAll(pDExt, speciei)
{
    pDExt.set
    (
        speciei,
        new scalarField(pD[speciei].primitiveField())
    );
    pDExt[speciei].setSize(nTotalCellsX);
}

if (nGhostCells > 0)
{
    // 值槽布局: rho, rhoU(3), rPsi, e, evk[nEv], pD[nPd]
    const label nEvX = evkExt.size();
    const label nPdX = pDExt.size();
    const label nValsX = 6 + nEvX + nPdX;

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // --- 发送: 按 ghostSendCells 打包本侧单元值 (即对端的 ghost 层) ---
    forAll(procPatchIDs, ppI)
    {
        const processorPolyPatch& procPatchX =
            refCast<const processorPolyPatch>
            (
                mesh.boundaryMesh()[procPatchIDs[ppI]]
            );

        const labelList& sendCellsX = ghostSendCells[ppI];
        scalarField sendBufX(sendCellsX.size()*nValsX);

        forAll(sendCellsX, scI)
        {
            const label c = sendCellsX[scI];
            label s = scI*nValsX;

            sendBufX[s++] = rhoI[c];
            sendBufX[s++] = rhoUI[c][0];
            sendBufX[s++] = rhoUI[c][1];
            sendBufX[s++] = rhoUI[c][2];
            sendBufX[s++] = rPsiI[c];
            sendBufX[s++] = eI[c];

            for (label k = 0; k < nEvX; ++k)
            {
                sendBufX[s++] = evkExt[k][c];
            }
            for (label k = 0; k < nPdX; ++k)
            {
                sendBufX[s++] = pDExt[k][c];
            }
        }

        UOPstream toNbr(procPatchX.neighbProcNo(), pBufs);
        toNbr << sendBufX;
    }

    pBufs.finishedSends();

    // --- 接收: 解包到扩展场 ghost 区 ---
    // 接收缓冲序 (face, d) 与本侧 ghost 编号一一对应
    forAll(procPatchIDs, ppI)
    {
        const processorPolyPatch& procPatchX =
            refCast<const processorPolyPatch>
            (
                mesh.boundaryMesh()[procPatchIDs[ppI]]
            );

        UIPstream fromNbr(procPatchX.neighbProcNo(), pBufs);
        scalarField recvBufX;
        fromNbr >> recvBufX;

        const label ghostBaseX = nLocalCellsX + ghostStart[ppI];
        const label nGhostsX = recvBufX.size()/nValsX;

        for (label gi = 0; gi < nGhostsX; ++gi)
        {
            const label g = ghostBaseX + gi;
            label s = gi*nValsX;

            rhoI[g] = recvBufX[s++];
            rhoUI[g] = vector(recvBufX[s], recvBufX[s+1], recvBufX[s+2]);
            s += 3;
            rPsiI[g] = recvBufX[s++];
            eI[g] = recvBufX[s++];

            for (label k = 0; k < nEvX; ++k)
            {
                evkExt[k][g] = recvBufX[s++];
            }
            for (label k = 0; k < nPdX; ++k)
            {
                pDExt[k][g] = recvBufX[s++];
            }
        }
    }
}

// ************************************************************************* //
