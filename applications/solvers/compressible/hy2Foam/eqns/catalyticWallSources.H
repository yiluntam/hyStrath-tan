//-----------------------------------------------------------
// Catalytic wall: boundary completion of the species-diffusive
// enthalpy flux into the energy equations
// 催化壁面：组分扩散焓通量的边界完成项（进入能量方程）
//-----------------------------------------------------------

// multiSpeciesHeatSource() / multiSpeciesVEHeatSource() fill the internal
// cells only; their boundary-face values are zero, so the enthalpy that
// the catalytic species fluxes transport through the wall faces never
// enters the energy equations. This include refreshes the boundary
// values of the completion fields constructed in createCatalyticFields.H
// so that, once added by the energy equations as
//     eqn += fvc::div(catalyticFlux*),   (same operator as the interior
//                                         species-enthalpy terms)
// the boundary contribution is present with a consistent discretisation.
//
// Conventions (identical to getDiffusiveWallHeatFlux()):
//   Is = -rhoD*snGrad(Y)  [kg/(m^2 s)], positive leaving the gas
//   Js = Is - Y_w*sum(Is) Stefan-corrected mass flux
// The gas side uses SENSIBLE enthalpies (hs, Hvels) mirroring the
// interior source terms (hyStrath's e is a sensible-energy formulation:
// the formation-enthalpy share of the recombination heat is not part of
// the gas energy budget). The wall REPORT (wallHeatFlux_cat) uses
// ABSOLUTE enthalpies (ha) so that the full recombination heat
// omega*Delta_h_recomb appears in the reported wall heat flux; the
// difference, omega*Delta_h_f, is carried by the gas-side chemical pool
// (composition change), consistently with the gas-phase reaction
// treatment.
//
// multiSpeciesHeatSource()/VEHeatSource() 只填充内部单元（边界面值
// 为零），催化组分通量穿过壁面所输运的焓因此从未进入能量方程。
// 本文件逐迭代刷新 createCatalyticFields.H 构造的完成通量场的边界
// 值。气侧用显焓（与内部源项一致，hyStrath 的 e 为显能形式），
// 壁面报告用绝对焓（复合反应热完整出现在报告热流中），差额
// ω·Δh_f 由气相化学池（组分变化）承担，与气相反应处理一致。

if (catalyticWallPresent)
{
    const labelList& heavyIds = composition.heavySpeciesIds();

    forAll(catalyticPatchIDs, i)
    {
        const label patchi = catalyticPatchIDs[i];
        const fvPatch& catPatch = mesh.boundary()[patchi];
        const label nFaces = catPatch.size();

        // Wall-face gas-side thermodynamic state (same fields as
        // getDiffusiveWallHeatFlux())
        const fvPatchScalarField& pp = p.boundaryField()[patchi];
        const fvPatchScalarField& pTt = thermo.T().boundaryField()[patchi];

        // Pass 1: raw boundary mass fluxes of the heavy species and
        // their sum (electrons excluded, catalysis acts on neutrals)
        // 第一遍：重组分的原始边界质量通量及其总和
        PtrList<scalarField> isList(Y.size());
        scalarField sumIs(nFaces, 0.0);

        forAll(heavyIds, hi)
        {
            const label speciei = heavyIds[hi];

            const fvPatchScalarField& prhoD =
                speciesDiffusion().rhoD(speciei)
                .boundaryField()[patchi];

            isList.set
            (
                speciei,
                new scalarField
                (
                    -prhoD*Y[speciei].boundaryField()[patchi].snGrad()
                )
            );

            sumIs += isList[speciei];
        }

        // Pass 2: Stefan correction and enthalpy-weighted sums
        // 第二遍：Stefan 修正与焓加权求和
        scalarField fFull(nFaces, 0.0);   // sum_s hs_s*Js_s,    sensible
        scalarField fVe(nFaces, 0.0);     // sum_s Hvels_s*Js_s, vibro-el.
        scalarField fCat(nFaces, 0.0);    // sum_s ha_s*Js_s,    wall report

        PtrList<scalarField> fVeK(noSolvedVibEq);
        for (label speciek = 0; speciek < noSolvedVibEq; speciek++)
        {
            fVeK.set(speciek, new scalarField(nFaces, 0.0));
        }

        forAll(heavyIds, hi)
        {
            const label speciei = heavyIds[hi];

            const fvPatchScalarField& pTv =
                composition.Tv()[speciei].boundaryField()[patchi];

            // Stefan correction, identical to the Y-equation treatment
            // 与 Y 方程一致的 Stefan 修正
            scalarField js = isList[speciei];
            if (!speciesDiffusion().useNonCorrected())
            {
                js -= Y[speciei].boundaryField()[patchi]*sumIs;
            }

            // Pointwise enthalpies at the wall faces
            // 壁面处的逐点焓值
            scalarField hsFace(nFaces, 0.0);
            scalarField hvelFace(nFaces, 0.0);
            scalarField haFace(nFaces, 0.0);

            forAll(js, facei)
            {
                hsFace[facei] =
                    speciesDiffusion().hs
                    (
                        speciei, pp[facei], pTt[facei], pTv[facei]
                    );
                hvelFace[facei] =
                    composition.Hvels
                    (
                        speciei, pp[facei], pTv[facei]
                    );
                haFace[facei] =
                    speciesDiffusion().ha
                    (
                        speciei, pp[facei], pTt[facei], pTv[facei]
                    );
            }

            fFull += hsFace*js;
            fVe += hvelFace*js;
            fCat += haFace*js;

            // Per-molecule-Tv share (same species indexing as the
            // per-species ev-equations, speciek < noSolvedVibEq)
            // 逐分子振动温度份额（与逐物种 ev 方程的物种索引一致）
            if (speciei < noSolvedVibEq)
            {
                fVeK[speciei] = hvelFace*js;
            }
        }

        // Assign the boundary values of the completion fields. The
        // vector value is (scalar flux)*(outward face normal), so that
        // fvc::div reproduces the missing boundary contribution with the
        // same operator sign as the interior species-enthalpy terms:
        // the gas loses the sensible enthalpy of the recombining atoms.
        // 赋完成通量场的边界面值：向量值 = 标量通量×外法向，使
        // fvc::div 以与内部源项一致的算符符号再现缺失的边界贡献
        // （气相失去复合原子的显焓）。
        vectorField fVec(nFaces, vector::zero);
        const vectorField nHat(catPatch.nf());

        if (downgradeSingleT)
        {
            forAll(fFull, facei)
            {
                fVec[facei] = fFull[facei]*nHat[facei];
            }
            catalyticFluxE->boundaryFieldRef()[patchi] = fVec;
        }
        else if (downgradeSingleTv)
        {
            forAll(fFull, facei)
            {
                fVec[facei] =
                    (fFull[facei] - fVe[facei])*nHat[facei];
            }
            catalyticFluxTr->boundaryFieldRef()[patchi] = fVec;

            forAll(fVe, facei)
            {
                fVec[facei] = fVe[facei]*nHat[facei];
            }
            catalyticFluxVe->boundaryFieldRef()[patchi] = fVec;
        }
        else
        {
            // One Tv per molecule: the et-share is the full sensible
            // share minus the vibro-electronic shares carried by the
            // per-species ev-equations
            // 逐分子 Tv：et 份额 = 全显焓份额 − 各 ev 方程承担的
            // 振电份额之和
            scalarField fVeSum(nFaces, 0.0);
            for (label speciek = 0; speciek < noSolvedVibEq; speciek++)
            {
                fVeSum += fVeK[speciek];
            }

            forAll(fFull, facei)
            {
                fVec[facei] =
                    (fFull[facei] - fVeSum[facei])*nHat[facei];
            }
            catalyticFluxTr->boundaryFieldRef()[patchi] = fVec;

            for (label speciek = 0; speciek < noSolvedVibEq; speciek++)
            {
                forAll(fVeK[speciek], facei)
                {
                    fVec[facei] = fVeK[speciek][facei]*nHat[facei];
                }
                catalyticFluxVeK->operator[](speciek)
                    .boundaryFieldRef()[patchi] = fVec;
            }
        }

        // Reported catalytic wall heat flux, positive into the wall
        // 报告用催化壁面热流（指向壁面为正）
        wallHeatFlux_cat.boundaryFieldRef()[patchi] = fCat;
    }
}
