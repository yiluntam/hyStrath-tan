<h1 align="center">hyStrath</h1>  

#### Hypersonic / Rarefied gas dynamics code developments under license GPL-3.0 
#### The only platform to conjointly host open-source CFD and DSMC codes designed for atmospheric re-entry analysis

#### The *Fleming* release includes  
+ *hyFoam*: a CFD solver for supersonic combusting flows   
+ *hy2Foam*: a CFD solver for hypersonic reacting flows with MHD capabilities  
+ *dsmcFoam+*: the direct simulation Monte Carlo (DSMC) code with all the latest features  
+ *pdFoam*: a hybrid PIC-DSMC solver   

#### Please visit the [_hyStrath_ website](https://hystrath.github.io/)

---

## hy2Foam 高阶格式扩展 (High-Order Schemes)

基于原始 [hyStrath](https://github.com/hystrath/hyStrath) 的 hy2Foam 求解器扩展, 添加以下功能:

### 通量格式与对流重构 (两级配置)

在 `system/fvSchemes` 中通过 `fluxScheme` 选择通量格式, `convSchemes` 选择
对流重构格式 (后者仅对 AUSM/Roe 生效):

| fluxScheme | 说明 |
|------------|------|
| `Kurganov` | 原生 KNP 中心格式 (默认), 不受新选项影响, 结果与原 hyStrath 一致 |
| `AUSM` | AUSM+-up 通量, 重构由 `convSchemes` 决定 |
| `Roe` | Roe 通量差分分裂 + Harten-Hyman 熵修正, 重构由 `convSchemes` 决定 |

| convSchemes | 精度 | 说明 |
|------------|------|------|
| `MUSCL` | 2–3 阶 | MUSCL 重构 (κ/限制器可调), 默认 |
| `WENOZ5` | 5 阶 | WENO-Z (Borges 2008) |
| `WENOZ7` | 7 阶 | WENO-Z |

配套 fvSchemes 可选字段 (仅在对应格式下生效, 其他格式指定了也不起作用):
- `MUSCLLimiter` (minmod | vanLeer | superbee | vanAlbada, 默认 vanAlbada) —
  仅 convSchemes = MUSCL 时生效的 TVD 限制器
- `MUSCLKappa` (标量, 默认 1/3) — 仅 convSchemes = MUSCL 时生效,
  1/3 → 3 阶, 0/±1 → 2 阶
- `EntropyCoeff` (标量, 默认 0.001) — 仅 fluxScheme = Roe 时生效,
  Harten-Hyman 熵修正阈值系数, delta = coeff·(|Un_rel|+c),
  |λ|<δ 时 |λ| ← (λ²+δ²)/(2δ)

旧关键字 `MUSCL_AUSM` / `WENOZ_AUSM` / `WENOZ7_AUSM` 仍被接受并自动映射
到 `fluxScheme AUSM + convSchemes MUSCL/WENOZ5/WENOZ7` (启动时打印提示).

注: AUSM/Roe 系列仅适用于结构化网格 (blockMesh 六面体), 非轴对齐面自动回退为 1 阶 upwind.

### 梯度格式

| gradScheme | 精度 | 说明 |
|------------|------|------|
| `sixth` | 6 阶 | 六阶紧致梯度 (需单独编译到 finiteVolume 库) |

### 新增/修改文件

```
hy2Foam/numerics/
├── AUSM_Flux.H                  # [新] AUSM+-up 通量函数
├── Roe_Flux.H                   # [新] Roe 通量差分分裂 (含熵修正)
├── MUSCL_Kernels.H              # [新] MUSCL 插值核 (κ/限制器可选)
├── WENOZ_Kernels.H              # [新] WENO-Z 插值核 (Borges 2008)
├── highOrderKernels.H           # [新] 模板/重构辅助函数
├── highOrderData.H              # [新] WENO 权重数据
├── MUSCL_WENO_Reconstruct.H     # [新] 每步面值重构
├── structuredStencilBuild.H     # [新] 结构化模板构建
├── readFluxScheme.H             # [改] 注册新通量方案, 读取 fvSchemes 选项
├── fluxesCalculation.H          # [改] 扩展通量条件
hy2Foam/
├── hy2Foam_include.H            # [改] 添加头文件引用
├── hy2Foam_solver.H             # [改] 插入重构步骤

OpenFOAM-v1706/.../gradSchemes/sixthGrad/
├── sixthGrad.H                  # [新] 六阶梯度类声明
├── sixthGrad.C                  # [新] 六阶梯度实现
└── sixthGrads.C                 # [新] 运行时注册

template/                        # [新] 算例模板
├── 0/                           # 边界条件模板
├── constant/                    # 物理模型模板
└── system/                      # 求解设置模板
```

### Bug 修复: 粘性能量方程的振动-平动交叉泄漏

**问题**: 原始 `eEqnViscous.H` 使用总热导率 `kappa = kappatr + kappave` 扩散总内能 `e = et + ev`,
然后 `et = e - ev` 提取平动能。此过程产生一个虚假的 `kappatr/Cv * grad(ev)` 热流项,
导致振动适应系数变化时平动热流出现非物理的 ~19% 差异。

**修复** (`eEqnViscous.H`): 直接对平动能 `et` 求解扩散方程, 使用平动热导率 `kappatr/CvtrMix`,
求解后通过 `e = et + ev` 重构总内能。振动能扩散完全由 `evEqnViscous` 处理 (`kappave`)。

**后续修复: 组分扩散焓源项的振动能双重计入**。改为求解 `et` 后, `eEqnViscous` 中的组分扩散
焓源项 `multiSpeciesHeatSource()` (输运全量焓 `hs = hts + hevel`) 未同步缩减——hevel 份额已由
`evEqnViscous` 的 `hevel·J` 项输运, 在总能预算中被计入两次。利用恒等式 `hs = hts + hevel`
(其中 `hvs == evs`、`hels == eels`, 因为焓的 pv 份额只属于平动模式), 修复按温度模式分派:
- `downgradeSingleT`: 保留全量 `sum_j hs_j * J_j` (无独立 ev 方程, 全量应归此处);
- `downgradeSingleTv`: 减去 `multiSpeciesVEHeatSource()` (与 ev 方程所加项严格互逆);
- 逐分子 Tv 模式: 减去 `sum_k hevel_k * J_k` (k 遍历 solvedVibEq 物种, 与 `evEqnViscous`
  的逐组分项逐项对应)。

修复后的算子分裂:
- `evEqnViscous`:  rho * d(ev)/dt = div(kappave * grad(Tv))  + div(hevel·J)
- `eEqnViscous`:   rho * d(et)/dt = div(kappatr * grad(Ttr)) + div(hts·J)
- 合计:           rho * d(e)/dt  = div(kappa * grad(T)) + div(hs·J)  ← 物理正确

### 原有功能完全保留

当 `fluxScheme` 为 `Kurganov` 或 `Tadmor` 时, 所有新代码无操作, 结果与原始 hyStrath 完全一致.
