//------------------------------------------------------------------------------
// MUSCL/WENO 高精度面重构 (每时间步) — 仅执行代码
//
// 在 upwindInterpolation.H 之后、fluxesCalculation.H 之前执行。
// 对轴对齐内部面覆写 rho_pos/neg, rhoU_pos/neg, rPsi_pos/neg, e_pos/neg,
// evk_pos/neg, pD_pos/neg, U_pos/neg, p_pos/neg, cSf_pos/neg 的内部面值。
// 边界面值保持不变 (来自 upwindInterpolation.H 的 BC 处理结果)。
//
// 辅助函数 (ghostStencilCell, stencilScalar, stencilVector, getDirLists)
// 定义在 numerics/highOrderKernels.H 中,
// MUSCL-3/WENO-Z 核函数定义在 numerics/MUSCL_Kernels.H, WENOZ_Kernels.H 中,
// 均由 hy2Foam_include.H 在文件作用域引入。
//
// 结构化模板数据 (axisFaces, cell_plus_x, etc.)
// 由 readFluxScheme.H 在函数作用域声明, structuredStencilBuild.H 初始化。
//------------------------------------------------------------------------------

// 仅在 MUSCL/WENO 方案激活时执行
if
(
    fluxScheme == "MUSCL_AUSM"
 or fluxScheme == "WENOZ_AUSM"
 or fluxScheme == "WENOZ7_AUSM"
)
{
    const surfaceScalarField& magSf = mesh.magSf();

    // 获取 vol 场内部值 (直接引用 primitiveField)
    const scalarField& rhoI   = rho.primitiveField();
    const vectorField& rhoUI  = rhoU.primitiveField();
    const scalarField& rPsiI  = rPsi.primitiveField();
    const scalarField& eI     = e.primitiveField();
    const scalarField& gammatrI = gammatr.primitiveField();

    // pos/neg 场的非 const 内部 field 引用
    scalarField& rho_posI   = rho_pos.primitiveFieldRef();
    scalarField& rho_negI   = rho_neg.primitiveFieldRef();
    vectorField& rhoU_posI  = rhoU_pos.primitiveFieldRef();
    vectorField& rhoU_negI  = rhoU_neg.primitiveFieldRef();
    scalarField& rPsi_posI  = rPsi_pos.primitiveFieldRef();
    scalarField& rPsi_negI  = rPsi_neg.primitiveFieldRef();
    scalarField& e_posI     = e_pos.primitiveFieldRef();
    scalarField& e_negI     = e_neg.primitiveFieldRef();
    vectorField& U_posI     = U_pos.primitiveFieldRef();
    vectorField& U_negI     = U_neg.primitiveFieldRef();
    scalarField& p_posI     = p_pos.primitiveFieldRef();
    scalarField& p_negI     = p_neg.primitiveFieldRef();
    scalarField& cSf_posI   = cSf_pos.primitiveFieldRef();
    scalarField& cSf_negI   = cSf_neg.primitiveFieldRef();

    // --- 循环所有内部面 ---
    forAll(mesh.owner(), facei)
    {
        const label dirTag = axisFaces[facei];
        if (dirTag == 0) continue;  // 非轴对齐面, 保留 upwind 值

        const label O = mesh.owner()[facei];
        const label N = mesh.neighbour()[facei];
        const scalar magSf_i = magSf[facei];

        // 获取方向邻居列表
        const labelList* plusDir  = nullptr;
        const labelList* minusDir = nullptr;
        getDirLists
        (
            dirTag,
            cell_plus_x, cell_minus_x,
            cell_plus_y, cell_minus_y,
            cell_plus_z, cell_minus_z,
            plusDir, minusDir
        );

        if (!plusDir || !minusDir) continue;

        // ============================================================
        // 标量重构: rho, rPsi, e
        // ============================================================
        scalar rhoL, rhoR, rPsiL, rPsiR, eL, eR;

        if (fluxScheme == "MUSCL_AUSM")
        {
            // MUSCL-3: {O-1, O, N=L+1, N+1=L+2}
            const scalar rho_m1 = stencilScalar(rhoI, O, -1, *plusDir, *minusDir);
            const scalar rho_p2 = stencilScalar(rhoI, N,  1, *plusDir, *minusDir);
            rhoL = muscl3FaceL(rho_m1, rhoI[O], rhoI[N]);
            rhoR = muscl3FaceR(rhoI[O], rhoI[N], rho_p2);

            const scalar rPsi_m1 = stencilScalar(rPsiI, O, -1, *plusDir, *minusDir);
            const scalar rPsi_p2 = stencilScalar(rPsiI, N,  1, *plusDir, *minusDir);
            rPsiL = muscl3FaceL(rPsi_m1, rPsiI[O], rPsiI[N]);
            rPsiR = muscl3FaceR(rPsiI[O], rPsiI[N], rPsi_p2);

            const scalar e_m1 = stencilScalar(eI, O, -1, *plusDir, *minusDir);
            const scalar e_p2 = stencilScalar(eI, N,  1, *plusDir, *minusDir);
            eL = muscl3FaceL(e_m1, eI[O], eI[N]);
            eR = muscl3FaceR(eI[O], eI[N], e_p2);
        }
        else if (fluxScheme == "WENOZ_AUSM")
        {
            // WENO-Z-5: L={O-2,O-1,O,N,N+1}, R={O-1,O,N,N+1,N+2}
            const scalar rho_m2 = stencilScalar(rhoI, O, -2, *plusDir, *minusDir);
            const scalar rho_m1 = stencilScalar(rhoI, O, -1, *plusDir, *minusDir);
            const scalar rho_p1 = stencilScalar(rhoI, N,  1, *plusDir, *minusDir);
            const scalar rho_p2 = stencilScalar(rhoI, N,  2, *plusDir, *minusDir);
            rhoL = wenoZ5_faceL(rho_m2, rho_m1, rhoI[O], rhoI[N], rho_p1);
            rhoR = wenoZ5_faceR(rho_m1, rhoI[O], rhoI[N], rho_p1, rho_p2);

            const scalar rPsi_m2 = stencilScalar(rPsiI, O, -2, *plusDir, *minusDir);
            const scalar rPsi_m1 = stencilScalar(rPsiI, O, -1, *plusDir, *minusDir);
            const scalar rPsi_p1 = stencilScalar(rPsiI, N,  1, *plusDir, *minusDir);
            const scalar rPsi_p2 = stencilScalar(rPsiI, N,  2, *plusDir, *minusDir);
            rPsiL = wenoZ5_faceL(rPsi_m2, rPsi_m1, rPsiI[O], rPsiI[N], rPsi_p1);
            rPsiR = wenoZ5_faceR(rPsi_m1, rPsiI[O], rPsiI[N], rPsi_p1, rPsi_p2);

            const scalar e_m2 = stencilScalar(eI, O, -2, *plusDir, *minusDir);
            const scalar e_m1 = stencilScalar(eI, O, -1, *plusDir, *minusDir);
            const scalar e_p1 = stencilScalar(eI, N,  1, *plusDir, *minusDir);
            const scalar e_p2 = stencilScalar(eI, N,  2, *plusDir, *minusDir);
            eL = wenoZ5_faceL(e_m2, e_m1, eI[O], eI[N], e_p1);
            eR = wenoZ5_faceR(e_m1, eI[O], eI[N], e_p1, e_p2);
        }
        else // WENOZ7_AUSM
        {
            // WENO-Z-7: L={O-3,O-2,O-1,O,N,N+1,N+2}, R={O-2,O-1,O,N,N+1,N+2,N+3}
            const scalar rho_m3 = stencilScalar(rhoI, O, -3, *plusDir, *minusDir);
            const scalar rho_m2 = stencilScalar(rhoI, O, -2, *plusDir, *minusDir);
            const scalar rho_m1 = stencilScalar(rhoI, O, -1, *plusDir, *minusDir);
            const scalar rho_p1 = stencilScalar(rhoI, N,  1, *plusDir, *minusDir);
            const scalar rho_p2 = stencilScalar(rhoI, N,  2, *plusDir, *minusDir);
            const scalar rho_p3 = stencilScalar(rhoI, N,  3, *plusDir, *minusDir);
            rhoL = wenoZ7_faceL(rho_m3, rho_m2, rho_m1, rhoI[O], rhoI[N], rho_p1, rho_p2);
            rhoR = wenoZ7_faceR(rho_m2, rho_m1, rhoI[O], rhoI[N], rho_p1, rho_p2, rho_p3);

            const scalar rPsi_m3 = stencilScalar(rPsiI, O, -3, *plusDir, *minusDir);
            const scalar rPsi_m2 = stencilScalar(rPsiI, O, -2, *plusDir, *minusDir);
            const scalar rPsi_m1 = stencilScalar(rPsiI, O, -1, *plusDir, *minusDir);
            const scalar rPsi_p1 = stencilScalar(rPsiI, N,  1, *plusDir, *minusDir);
            const scalar rPsi_p2 = stencilScalar(rPsiI, N,  2, *plusDir, *minusDir);
            const scalar rPsi_p3 = stencilScalar(rPsiI, N,  3, *plusDir, *minusDir);
            rPsiL = wenoZ7_faceL(rPsi_m3, rPsi_m2, rPsi_m1, rPsiI[O], rPsiI[N], rPsi_p1, rPsi_p2);
            rPsiR = wenoZ7_faceR(rPsi_m2, rPsi_m1, rPsiI[O], rPsiI[N], rPsi_p1, rPsi_p2, rPsi_p3);

            const scalar e_m3 = stencilScalar(eI, O, -3, *plusDir, *minusDir);
            const scalar e_m2 = stencilScalar(eI, O, -2, *plusDir, *minusDir);
            const scalar e_m1 = stencilScalar(eI, O, -1, *plusDir, *minusDir);
            const scalar e_p1 = stencilScalar(eI, N,  1, *plusDir, *minusDir);
            const scalar e_p2 = stencilScalar(eI, N,  2, *plusDir, *minusDir);
            const scalar e_p3 = stencilScalar(eI, N,  3, *plusDir, *minusDir);
            eL = wenoZ7_faceL(e_m3, e_m2, e_m1, eI[O], eI[N], e_p1, e_p2);
            eR = wenoZ7_faceR(e_m2, e_m1, eI[O], eI[N], e_p1, e_p2, e_p3);
        }

        // 正性保持: 面值必须为正, 否则回退到 upwind (cell 中心值)
        if (rhoL <= 0)  rhoL  = rhoI[O];
        if (rhoR <= 0)  rhoR  = rhoI[N];
        if (rPsiL <= 0) rPsiL = rPsiI[O];
        if (rPsiR <= 0) rPsiR = rPsiI[N];
        if (eL <= 0)    eL    = eI[O];
        if (eR <= 0)    eR    = eI[N];

        // 覆写标量面场
        rho_posI[facei]  = rhoL;
        rho_negI[facei]  = rhoR;
        rPsi_posI[facei] = rPsiL;
        rPsi_negI[facei] = rPsiR;
        e_posI[facei]    = eL;
        e_negI[facei]    = eR;

        // ============================================================
        // 矢量重构: rhoU (动量)
        // ============================================================
        vector rhoUL, rhoUR;

        if (fluxScheme == "MUSCL_AUSM")
        {
            const vector rhoU_m1 = stencilVector(rhoUI, O, -1, *plusDir, *minusDir);
            const vector rhoU_p2 = stencilVector(rhoUI, N,  1, *plusDir, *minusDir);
            rhoUL = muscl3FaceL(rhoU_m1, rhoUI[O], rhoUI[N]);
            rhoUR = muscl3FaceR(rhoUI[O], rhoUI[N], rhoU_p2);
        }
        else if (fluxScheme == "WENOZ_AUSM")
        {
            const vector rhoU_m2 = stencilVector(rhoUI, O, -2, *plusDir, *minusDir);
            const vector rhoU_m1 = stencilVector(rhoUI, O, -1, *plusDir, *minusDir);
            const vector rhoU_p1 = stencilVector(rhoUI, N,  1, *plusDir, *minusDir);
            const vector rhoU_p2 = stencilVector(rhoUI, N,  2, *plusDir, *minusDir);

            // WENO 逐分量应用
            for (label cmp = 0; cmp < 3; ++cmp)
            {
                rhoUL[cmp] = wenoZ5_faceL
                (
                    rhoU_m2[cmp], rhoU_m1[cmp], rhoUI[O][cmp],
                    rhoUI[N][cmp], rhoU_p1[cmp]
                );
                rhoUR[cmp] = wenoZ5_faceR
                (
                    rhoU_m1[cmp], rhoUI[O][cmp], rhoUI[N][cmp],
                    rhoU_p1[cmp], rhoU_p2[cmp]
                );
            }
        }
        else // WENOZ7_AUSM
        {
            const vector rhoU_m3 = stencilVector(rhoUI, O, -3, *plusDir, *minusDir);
            const vector rhoU_m2 = stencilVector(rhoUI, O, -2, *plusDir, *minusDir);
            const vector rhoU_m1 = stencilVector(rhoUI, O, -1, *plusDir, *minusDir);
            const vector rhoU_p1 = stencilVector(rhoUI, N,  1, *plusDir, *minusDir);
            const vector rhoU_p2 = stencilVector(rhoUI, N,  2, *plusDir, *minusDir);
            const vector rhoU_p3 = stencilVector(rhoUI, N,  3, *plusDir, *minusDir);

            for (label cmp = 0; cmp < 3; ++cmp)
            {
                rhoUL[cmp] = wenoZ7_faceL
                (
                    rhoU_m3[cmp], rhoU_m2[cmp], rhoU_m1[cmp], rhoUI[O][cmp],
                    rhoUI[N][cmp], rhoU_p1[cmp], rhoU_p2[cmp]
                );
                rhoUR[cmp] = wenoZ7_faceR
                (
                    rhoU_m2[cmp], rhoU_m1[cmp], rhoUI[O][cmp], rhoUI[N][cmp],
                    rhoU_p1[cmp], rhoU_p2[cmp], rhoU_p3[cmp]
                );
            }
        }

        rhoU_posI[facei] = rhoUL;
        rhoU_negI[facei] = rhoUR;

        // ============================================================
        // 振动能量重构
        // ============================================================
        if (!downgradeSingleT)
        {
            if (downgradeSingleTv)
            {
                // 单振动模式: evk_pos[0] / evk_neg[0]
                const scalarField& evI = ev.primitiveField();
                scalarField& evk_pos0 = evk_pos[0].primitiveFieldRef();
                scalarField& evk_neg0 = evk_neg[0].primitiveFieldRef();

                if (fluxScheme == "MUSCL_AUSM")
                {
                    const scalar ev_m1 = stencilScalar(evI, O, -1, *plusDir, *minusDir);
                    const scalar ev_p2 = stencilScalar(evI, N,  1, *plusDir, *minusDir);
                    evk_pos0[facei] = muscl3FaceL(ev_m1, evI[O], evI[N]);
                    evk_neg0[facei] = muscl3FaceR(evI[O], evI[N], ev_p2);
                }
                else if (fluxScheme == "WENOZ_AUSM")
                {
                    const scalar ev_m2 = stencilScalar(evI, O, -2, *plusDir, *minusDir);
                    const scalar ev_m1 = stencilScalar(evI, O, -1, *plusDir, *minusDir);
                    const scalar ev_p1 = stencilScalar(evI, N,  1, *plusDir, *minusDir);
                    const scalar ev_p2 = stencilScalar(evI, N,  2, *plusDir, *minusDir);
                    evk_pos0[facei] = wenoZ5_faceL(ev_m2, ev_m1, evI[O], evI[N], ev_p1);
                    evk_neg0[facei] = wenoZ5_faceR(ev_m1, evI[O], evI[N], ev_p1, ev_p2);
                }
                else
                {
                    const scalar ev_m3 = stencilScalar(evI, O, -3, *plusDir, *minusDir);
                    const scalar ev_m2 = stencilScalar(evI, O, -2, *plusDir, *minusDir);
                    const scalar ev_m1 = stencilScalar(evI, O, -1, *plusDir, *minusDir);
                    const scalar ev_p1 = stencilScalar(evI, N,  1, *plusDir, *minusDir);
                    const scalar ev_p2 = stencilScalar(evI, N,  2, *plusDir, *minusDir);
                    const scalar ev_p3 = stencilScalar(evI, N,  3, *plusDir, *minusDir);
                    evk_pos0[facei] = wenoZ7_faceL(ev_m3, ev_m2, ev_m1, evI[O], evI[N], ev_p1, ev_p2);
                    evk_neg0[facei] = wenoZ7_faceR(ev_m2, ev_m1, evI[O], evI[N], ev_p1, ev_p2, ev_p3);
                }
            }
            else
            {
                // 多物种: evk_pos[speciei] / evk_neg[speciei]
                forAll(evk_pos, speciei)
                {
                    const scalarField& evkI = hevk[speciei].primitiveField();
                    scalarField& evk_posI_sp = evk_pos[speciei].primitiveFieldRef();
                    scalarField& evk_negI_sp = evk_neg[speciei].primitiveFieldRef();

                    if (fluxScheme == "MUSCL_AUSM")
                    {
                        const scalar ev_m1 = stencilScalar(evkI, O, -1, *plusDir, *minusDir);
                        const scalar ev_p2 = stencilScalar(evkI, N,  1, *plusDir, *minusDir);
                        evk_posI_sp[facei] = muscl3FaceL(ev_m1, evkI[O], evkI[N]);
                        evk_negI_sp[facei] = muscl3FaceR(evkI[O], evkI[N], ev_p2);
                    }
                    else if (fluxScheme == "WENOZ_AUSM")
                    {
                        const scalar ev_m2 = stencilScalar(evkI, O, -2, *plusDir, *minusDir);
                        const scalar ev_m1 = stencilScalar(evkI, O, -1, *plusDir, *minusDir);
                        const scalar ev_p1 = stencilScalar(evkI, N,  1, *plusDir, *minusDir);
                        const scalar ev_p2 = stencilScalar(evkI, N,  2, *plusDir, *minusDir);
                        evk_posI_sp[facei] = wenoZ5_faceL(ev_m2, ev_m1, evkI[O], evkI[N], ev_p1);
                        evk_negI_sp[facei] = wenoZ5_faceR(ev_m1, evkI[O], evkI[N], ev_p1, ev_p2);
                    }
                    else
                    {
                        const scalar ev_m3 = stencilScalar(evkI, O, -3, *plusDir, *minusDir);
                        const scalar ev_m2 = stencilScalar(evkI, O, -2, *plusDir, *minusDir);
                        const scalar ev_m1 = stencilScalar(evkI, O, -1, *plusDir, *minusDir);
                        const scalar ev_p1 = stencilScalar(evkI, N,  1, *plusDir, *minusDir);
                        const scalar ev_p2 = stencilScalar(evkI, N,  2, *plusDir, *minusDir);
                        const scalar ev_p3 = stencilScalar(evkI, N,  3, *plusDir, *minusDir);
                        evk_posI_sp[facei] = wenoZ7_faceL(ev_m3, ev_m2, ev_m1, evkI[O], evkI[N], ev_p1, ev_p2);
                        evk_negI_sp[facei] = wenoZ7_faceR(ev_m2, ev_m1, evkI[O], evkI[N], ev_p1, ev_p2, ev_p3);
                    }
                }

                // 部分密度 pD 重构
                forAll(pD_pos, speciei)
                {
                    const scalarField& pDI = pD[speciei].primitiveField();
                    scalarField& pD_posI_sp = pD_pos[speciei].primitiveFieldRef();
                    scalarField& pD_negI_sp = pD_neg[speciei].primitiveFieldRef();

                    if (fluxScheme == "MUSCL_AUSM")
                    {
                        const scalar pd_m1 = stencilScalar(pDI, O, -1, *plusDir, *minusDir);
                        const scalar pd_p2 = stencilScalar(pDI, N,  1, *plusDir, *minusDir);
                        pD_posI_sp[facei] = muscl3FaceL(pd_m1, pDI[O], pDI[N]);
                        pD_negI_sp[facei] = muscl3FaceR(pDI[O], pDI[N], pd_p2);
                    }
                    else if (fluxScheme == "WENOZ_AUSM")
                    {
                        const scalar pd_m2 = stencilScalar(pDI, O, -2, *plusDir, *minusDir);
                        const scalar pd_m1 = stencilScalar(pDI, O, -1, *plusDir, *minusDir);
                        const scalar pd_p1 = stencilScalar(pDI, N,  1, *plusDir, *minusDir);
                        const scalar pd_p2 = stencilScalar(pDI, N,  2, *plusDir, *minusDir);
                        pD_posI_sp[facei] = wenoZ5_faceL(pd_m2, pd_m1, pDI[O], pDI[N], pd_p1);
                        pD_negI_sp[facei] = wenoZ5_faceR(pd_m1, pDI[O], pDI[N], pd_p1, pd_p2);
                    }
                    else
                    {
                        const scalar pd_m3 = stencilScalar(pDI, O, -3, *plusDir, *minusDir);
                        const scalar pd_m2 = stencilScalar(pDI, O, -2, *plusDir, *minusDir);
                        const scalar pd_m1 = stencilScalar(pDI, O, -1, *plusDir, *minusDir);
                        const scalar pd_p1 = stencilScalar(pDI, N,  1, *plusDir, *minusDir);
                        const scalar pd_p2 = stencilScalar(pDI, N,  2, *plusDir, *minusDir);
                        const scalar pd_p3 = stencilScalar(pDI, N,  3, *plusDir, *minusDir);
                        pD_posI_sp[facei] = wenoZ7_faceL(pd_m3, pd_m2, pd_m1, pDI[O], pDI[N], pd_p1, pd_p2);
                        pD_negI_sp[facei] = wenoZ7_faceR(pd_m2, pd_m1, pDI[O], pDI[N], pd_p1, pd_p2, pd_p3);
                    }
                }
            }
        }

        // ============================================================
        // 派生场: U = rhoU / rho, p = rho * rPsi, cSf = sqrt(gamma*rPsi)*|Sf|
        // ============================================================
        U_posI[facei] = rhoUL / (rhoL + SMALL);
        U_negI[facei] = rhoUR / (rhoR + SMALL);

        p_posI[facei] = rhoL * rPsiL;
        p_negI[facei] = rhoR * rPsiR;

        cSf_posI[facei] = ::sqrt(gammatrI[O] * rPsiL) * magSf_i;
        cSf_negI[facei] = ::sqrt(gammatrI[N] * rPsiR) * magSf_i;
    }

    Info<< "MUSCL/WENO reconstruction applied on " << nAxisFaces
        << " axis-aligned faces." << endl;
}

// ************************************************************************* //
