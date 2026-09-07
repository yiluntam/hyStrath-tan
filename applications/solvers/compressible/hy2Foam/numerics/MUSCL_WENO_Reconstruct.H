//------------------------------------------------------------------------------
// MUSCL/WENO 高精度面重构 (每时间步) — 仅执行代码
//
// 在 upwindInterpolation.H 之后、fluxesCalculation.H 之前执行。
// 对轴对齐内部面覆写 rho_pos/neg, rhoU_pos/neg, rPsi_pos/neg, e_pos/neg,
// evk_pos/neg, pD_pos/neg, U_pos/neg, p_pos/neg, cSf_pos/neg 的内部面值。
// 边界面值保持不变 (来自 upwindInterpolation.H 的 BC 处理结果)。
//
// bad_recon 整体回退 (SU2 CEulerSolver::Upwind_Residual):
//   每面重构完成后检查两侧 ρ/p/e 与 Roe 平均声速²; 任一非物理
//   (或 NaN) 则该面 [全部变量] 两侧一起回退为单元中心值 (局部一阶),
//   而非旧版的逐变量修补。正常计算仅启动瞬态触发。
//
// 模板取值使用 "本地 + ghost" 扩展场 (ghostFieldExchange.H 组装并完成
// 处理器间多层交换), 使高阶模板在分区边界取得远端真实单元值。
// 扩展后的方向链在处理器边界处延伸到 ghost 单元
// (ghostStencilBuild.H), 走链代码本身无需感知并行。
//
// 辅助函数 (ghostStencilCell, stencilScalar, stencilVector, getDirLists)
// 定义在 numerics/highOrderKernels.H 中,
// MUSCL-3/WENO-Z 核函数定义在 numerics/MUSCL_Kernels.H, WENOZ_Kernels.H 中,
// 均由 hy2Foam_include.H 在文件作用域引入。
//
// 结构化模板数据 (axisFaces, cell_±_x/y/z, ghost 层映射)
// 由 readFluxScheme.H 在函数作用域声明,
// structuredStencilBuild.H + ghostStencilBuild.H 初始化。
//------------------------------------------------------------------------------

// 仅在 fluxScheme = AUSM/Roe 时执行 (重构方式由 convSchemes 选择)
if (fluxScheme == "AUSM" or fluxScheme == "Roe")
{
    const surfaceScalarField& magSf = mesh.magSf();

    // 获取 vol 场内部值 (本地 + ghost 扩展场; 串行时即原场拷贝)
    #include "ghostFieldExchange.H"

    // gammatr 仅在面两侧本地单元 O/N 取值, 无需 ghost 扩展
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

    // --- 重构暂存 (循环外分配一次; bad_recon 整体回退需要先算后写) ---
    scalar evL0 = 0.0;              // 单振动模式: 面左/右比振动能
    scalar evR0 = 0.0;
    scalarList evkL, evkR;          // 多物种: 比振动能
    scalarList pDL, pDR;            // 多物种: 部分密度
    if (!downgradeSingleT)
    {
        if (!downgradeSingleTv)
        {
            evkL.resize(evkExt.size());
            evkR.resize(evkExt.size());
            pDL.resize(pDExt.size());
            pDR.resize(pDExt.size());
        }
    }
    label badReconFaces = 0;

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

        if (convSchemes == "MUSCL")
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
        else if (convSchemes == "WENOZ5")
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
        else // WENOZ7
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

        // ============================================================
        // 矢量重构: rhoU (动量)
        // ============================================================
        vector rhoUL, rhoUR;

        if (convSchemes == "MUSCL")
        {
            const vector rhoU_m1 = stencilVector(rhoUI, O, -1, *plusDir, *minusDir);
            const vector rhoU_p2 = stencilVector(rhoUI, N,  1, *plusDir, *minusDir);
            rhoUL = muscl3FaceL(rhoU_m1, rhoUI[O], rhoUI[N]);
            rhoUR = muscl3FaceR(rhoUI[O], rhoUI[N], rhoU_p2);
        }
        else if (convSchemes == "WENOZ5")
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
        else // WENOZ7
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

        // ============================================================
        // 振动能量重构 (先算入暂存, 统一写入在 bad_recon 检查之后)
        // ============================================================
        if (!downgradeSingleT)
        {
            if (downgradeSingleTv)
            {
                // 单振动模式: evk_pos[0] / evk_neg[0]
                // (evI 为 ghostFieldExchange.H 的扩展场 evkExt[0])
                const scalarField& evI = evkExt[0];

                if (convSchemes == "MUSCL")
                {
                    const scalar ev_m1 = stencilScalar(evI, O, -1, *plusDir, *minusDir);
                    const scalar ev_p2 = stencilScalar(evI, N,  1, *plusDir, *minusDir);
                    evL0 = muscl3FaceL(ev_m1, evI[O], evI[N]);
                    evR0 = muscl3FaceR(evI[O], evI[N], ev_p2);
                }
                else if (convSchemes == "WENOZ5")
                {
                    const scalar ev_m2 = stencilScalar(evI, O, -2, *plusDir, *minusDir);
                    const scalar ev_m1 = stencilScalar(evI, O, -1, *plusDir, *minusDir);
                    const scalar ev_p1 = stencilScalar(evI, N,  1, *plusDir, *minusDir);
                    const scalar ev_p2 = stencilScalar(evI, N,  2, *plusDir, *minusDir);
                    evL0 = wenoZ5_faceL(ev_m2, ev_m1, evI[O], evI[N], ev_p1);
                    evR0 = wenoZ5_faceR(ev_m1, evI[O], evI[N], ev_p1, ev_p2);
                }
                else
                {
                    const scalar ev_m3 = stencilScalar(evI, O, -3, *plusDir, *minusDir);
                    const scalar ev_m2 = stencilScalar(evI, O, -2, *plusDir, *minusDir);
                    const scalar ev_m1 = stencilScalar(evI, O, -1, *plusDir, *minusDir);
                    const scalar ev_p1 = stencilScalar(evI, N,  1, *plusDir, *minusDir);
                    const scalar ev_p2 = stencilScalar(evI, N,  2, *plusDir, *minusDir);
                    const scalar ev_p3 = stencilScalar(evI, N,  3, *plusDir, *minusDir);
                    evL0 = wenoZ7_faceL(ev_m3, ev_m2, ev_m1, evI[O], evI[N], ev_p1, ev_p2);
                    evR0 = wenoZ7_faceR(ev_m2, ev_m1, evI[O], evI[N], ev_p1, ev_p2, ev_p3);
                }
            }
            else
            {
                // 多物种: evk / pD 逐组重构
                // (evkI/pDI 为 ghostFieldExchange.H 的扩展场 evkExt/pDExt)
                forAll(evkExt, speciei)
                {
                    const scalarField& evkI = evkExt[speciei];

                    if (convSchemes == "MUSCL")
                    {
                        const scalar ev_m1 = stencilScalar(evkI, O, -1, *plusDir, *minusDir);
                        const scalar ev_p2 = stencilScalar(evkI, N,  1, *plusDir, *minusDir);
                        evkL[speciei] = muscl3FaceL(ev_m1, evkI[O], evkI[N]);
                        evkR[speciei] = muscl3FaceR(evkI[O], evkI[N], ev_p2);
                    }
                    else if (convSchemes == "WENOZ5")
                    {
                        const scalar ev_m2 = stencilScalar(evkI, O, -2, *plusDir, *minusDir);
                        const scalar ev_m1 = stencilScalar(evkI, O, -1, *plusDir, *minusDir);
                        const scalar ev_p1 = stencilScalar(evkI, N,  1, *plusDir, *minusDir);
                        const scalar ev_p2 = stencilScalar(evkI, N,  2, *plusDir, *minusDir);
                        evkL[speciei] = wenoZ5_faceL(ev_m2, ev_m1, evkI[O], evkI[N], ev_p1);
                        evkR[speciei] = wenoZ5_faceR(ev_m1, evkI[O], evkI[N], ev_p1, ev_p2);
                    }
                    else
                    {
                        const scalar ev_m3 = stencilScalar(evkI, O, -3, *plusDir, *minusDir);
                        const scalar ev_m2 = stencilScalar(evkI, O, -2, *plusDir, *minusDir);
                        const scalar ev_m1 = stencilScalar(evkI, O, -1, *plusDir, *minusDir);
                        const scalar ev_p1 = stencilScalar(evkI, N,  1, *plusDir, *minusDir);
                        const scalar ev_p2 = stencilScalar(evkI, N,  2, *plusDir, *minusDir);
                        const scalar ev_p3 = stencilScalar(evkI, N,  3, *plusDir, *minusDir);
                        evkL[speciei] = wenoZ7_faceL(ev_m3, ev_m2, ev_m1, evkI[O], evkI[N], ev_p1, ev_p2);
                        evkR[speciei] = wenoZ7_faceR(ev_m2, ev_m1, evkI[O], evkI[N], ev_p1, ev_p2, ev_p3);
                    }
                }

                // 部分密度 pD 重构
                forAll(pDExt, speciei)
                {
                    const scalarField& pDI = pDExt[speciei];

                    if (convSchemes == "MUSCL")
                    {
                        const scalar pd_m1 = stencilScalar(pDI, O, -1, *plusDir, *minusDir);
                        const scalar pd_p2 = stencilScalar(pDI, N,  1, *plusDir, *minusDir);
                        pDL[speciei] = muscl3FaceL(pd_m1, pDI[O], pDI[N]);
                        pDR[speciei] = muscl3FaceR(pDI[O], pDI[N], pd_p2);
                    }
                    else if (convSchemes == "WENOZ5")
                    {
                        const scalar pd_m2 = stencilScalar(pDI, O, -2, *plusDir, *minusDir);
                        const scalar pd_m1 = stencilScalar(pDI, O, -1, *plusDir, *minusDir);
                        const scalar pd_p1 = stencilScalar(pDI, N,  1, *plusDir, *minusDir);
                        const scalar pd_p2 = stencilScalar(pDI, N,  2, *plusDir, *minusDir);
                        pDL[speciei] = wenoZ5_faceL(pd_m2, pd_m1, pDI[O], pDI[N], pd_p1);
                        pDR[speciei] = wenoZ5_faceR(pd_m1, pDI[O], pDI[N], pd_p1, pd_p2);
                    }
                    else
                    {
                        const scalar pd_m3 = stencilScalar(pDI, O, -3, *plusDir, *minusDir);
                        const scalar pd_m2 = stencilScalar(pDI, O, -2, *plusDir, *minusDir);
                        const scalar pd_m1 = stencilScalar(pDI, O, -1, *plusDir, *minusDir);
                        const scalar pd_p1 = stencilScalar(pDI, N,  1, *plusDir, *minusDir);
                        const scalar pd_p2 = stencilScalar(pDI, N,  2, *plusDir, *minusDir);
                        const scalar pd_p3 = stencilScalar(pDI, N,  3, *plusDir, *minusDir);
                        pDL[speciei] = wenoZ7_faceL(pd_m3, pd_m2, pd_m1, pDI[O], pDI[N], pd_p1, pd_p2);
                        pDR[speciei] = wenoZ7_faceR(pd_m2, pd_m1, pDI[O], pDI[N], pd_p1, pd_p2, pd_p3);
                    }
                }
            }
        }

        // ============================================================
        // bad_recon 整体检查 (SU2 CEulerSolver::Upwind_Residual):
        //   任一侧 ρ/p (=ρ·rPsi)/e 非物理, 或两侧重构值的 Roe 平均
        //   声速² < 0 → 该面所有变量两侧一起回退单元中心值 (一阶)。
        //   !(x > 0) 写法同时捕获 NaN。
        // ============================================================
        bool badRecon = !(rhoL > 0.0) || !(rhoR > 0.0)
                     || !(rPsiL > 0.0) || !(rPsiR > 0.0)
                     || !(eL > 0.0) || !(eR > 0.0);

        if (!badRecon)
        {
            // Roe 平均声速检查 (两温: 与 Roe_Flux.H 的 ĉ² 公式一致;
            // 侧值 γ = ρc²/p = gammatr, 因 c² = γtr·rPsi 且 p = ρ·rPsi)
            const vector UL = rhoUL/(rhoL + SMALL);
            const vector UR = rhoUR/(rhoR + SMALL);
            const scalar pL = rhoL*rPsiL;
            const scalar pR = rhoR*rPsiR;
            const scalar HL = eL + pL/rhoL + 0.5*magSqr(UL);
            const scalar HR = eR + pR/rhoR + 0.5*magSqr(UR);

            scalar evL = 0.0;
            scalar evR = 0.0;
            if (!downgradeSingleT)
            {
                if (downgradeSingleTv)
                {
                    evL = evL0;
                    evR = evR0;
                }
                else
                {
                    forAll(evkL, speciei)
                    {
                        evL += (pDL[speciei]/rhoL)*evkL[speciei];
                        evR += (pDR[speciei]/rhoR)*evkR[speciei];
                    }
                }
            }

            const scalar sL = Foam::sqrt(rhoL);
            const scalar sR = Foam::sqrt(rhoR);
            const scalar wL = sL/(sL + sR);
            const scalar wR = 1.0 - wL;

            const vector UHatC = wL*UL + wR*UR;
            const scalar HHatC = wL*HL + wR*HR;
            const scalar evHatC = wL*evL + wR*evR;
            const scalar gammaHatC = wL*gammatrI[O] + wR*gammatrI[N];

            const scalar cHat2C =
                (gammaHatC - 1.0)*(HHatC - 0.5*magSqr(UHatC) - evHatC);

            if (!(cHat2C >= 0.0))
            {
                badRecon = true;
            }
        }

        if (badRecon)
        {
            ++badReconFaces;
            rhoL = rhoI[O];      rhoR = rhoI[N];
            rPsiL = rPsiI[O];    rPsiR = rPsiI[N];
            eL = eI[O];          eR = eI[N];
            rhoUL = rhoUI[O];    rhoUR = rhoUI[N];

            if (!downgradeSingleT)
            {
                if (downgradeSingleTv)
                {
                    evL0 = evkExt[0][O];
                    evR0 = evkExt[0][N];
                }
                else
                {
                    forAll(evkL, speciei)
                    {
                        evkL[speciei] = evkExt[speciei][O];
                        evkR[speciei] = evkExt[speciei][N];
                    }
                    forAll(pDL, speciei)
                    {
                        pDL[speciei] = pDExt[speciei][O];
                        pDR[speciei] = pDExt[speciei][N];
                    }
                }
            }
        }

        // ============================================================
        // 统一写入面场
        // ============================================================
        rho_posI[facei]  = rhoL;
        rho_negI[facei]  = rhoR;
        rPsi_posI[facei] = rPsiL;
        rPsi_negI[facei] = rPsiR;
        e_posI[facei]    = eL;
        e_negI[facei]    = eR;

        rhoU_posI[facei] = rhoUL;
        rhoU_negI[facei] = rhoUR;

        if (!downgradeSingleT)
        {
            if (downgradeSingleTv)
            {
                scalarField& evk_pos0 = evk_pos[0].primitiveFieldRef();
                scalarField& evk_neg0 = evk_neg[0].primitiveFieldRef();
                evk_pos0[facei] = evL0;
                evk_neg0[facei] = evR0;
            }
            else
            {
                forAll(evkL, speciei)
                {
                    scalarField& evk_posI_sp = evk_pos[speciei].primitiveFieldRef();
                    scalarField& evk_negI_sp = evk_neg[speciei].primitiveFieldRef();
                    evk_posI_sp[facei] = evkL[speciei];
                    evk_negI_sp[facei] = evkR[speciei];
                }
                forAll(pDL, speciei)
                {
                    scalarField& pD_posI_sp = pD_pos[speciei].primitiveFieldRef();
                    scalarField& pD_negI_sp = pD_neg[speciei].primitiveFieldRef();
                    pD_posI_sp[facei] = pDL[speciei];
                    pD_negI_sp[facei] = pDR[speciei];
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

    if (badReconFaces > 0)
    {
        Info<< "MUSCL/WENO bad_recon faces (reverted to 1st order): "
            << badReconFaces << endl;
    }
}

// ************************************************************************* //
