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
- `Mach_inf` (标量, 必须为正) — 仅 fluxScheme = AUSM 时必需 (缺省启动报错),
  AUSM+-up 低马赫修正参数 fa/Mref 的远场马赫数, 取算例来流马赫数
- `EntropyCoeff` (标量, 默认 1.0) — 仅 fluxScheme = Roe 时生效,
  Harten-Hyman (1983) 逐波熵修正带宽系数:
  ε_k = EntropyCoeff·4·max(0, max(λ̂_k − λ_k,L, λ_k,R − λ̂_k)),
  |λ̂_k| < ε_k 时 |λ̂_k| ← (λ̂_k² + ε_k²)/(2ε_k)。
  压缩波/接触间断处 ε_k = 0, 修正自动关闭。默认 1.0 = SU2 NEMO 原样
  (含其固有的 4 倍安全因子); 取 0.25 可去掉该因子回到教科书 Harten-Hyman

旧关键字 `MUSCL_AUSM` / `WENOZ_AUSM` / `WENOZ7_AUSM` 仍被接受并自动映射
到 `fluxScheme AUSM + convSchemes MUSCL/WENOZ5/WENOZ7` (启动时打印提示).

注: AUSM/Roe 系列仅适用于结构化网格 (blockMesh 六面体); 存在非轴对齐
内部面时启动即报错退出 (structuredStencilBuild 要求 100% 轴对齐,
不再静默降阶). 非结构网格请用 Kurganov.

**并行运行**: MUSCL/WENO 重构已支持处理器边界 ghost 层交换
(`ghostStencilBuild.H` 初始化 + `ghostFieldExchange.H` 每步交换) —
方向链在分区边界延伸到远端进程单元 (深度 MUSCL=1 / WENOZ5=2 / WENOZ7=3),
高阶模板在分区边界不再退化为常值外推。串行结果不受影响;
`fluxScheme Kurganov/Tadmor` 路径完全不经过该机制。
限制: cyclic/periodic 边界与动网格暂不支持 (仍为常值外推/初始化期建链)。

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
├── ghostStencilBuild.H          # [新] 处理器边界 ghost 链构建 (初始化)
├── ghostFieldExchange.H         # [新] ghost 层场值交换 (每步)
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

### 已知限制 (2026-09 数值格式审查)

- **动网格路径无高阶重构**: `DyM/hy2DyMFoam_solver.H` 未 include
  `MUSCL_WENO_Reconstruct.H`, 动网格 + AUSM/Roe 时面状态静默保持 1 阶
  迎风 (通量函数的 ALE 分支正常, 但重构不执行, 亦无提示)。
- **重构状态无 pD/evk 正性检查**: bad_recon 仅检查 ρ/p/e 与 Roe 平均声速²,
  强激波附近组分部分密度 pD 与振动能 evk 重构可能为负并进入振动能通量
  (SU2 重构守恒变量, 无此暴露面)。
- **LTS 为半成品**: 反应源/温度时标限制在 `LTS/setrDeltaT.H` 中被禁用
  (作者标注 UNSTABLE), controlDict 的 `rDeltaTSmoothingCoeff` 为死参数
  (代码固定 0.1), `rDeltaTDampingCoeff`/`alphaTemp` 未使用 —
  反应流算例勿依赖 LTS 的化学刚性保护。

### 原有功能完全保留

当 `fluxScheme` 为 `Kurganov` 或 `Tadmor` 时, 所有新代码无操作, 结果与原始 hyStrath 完全一致.

## 催化壁面边界条件 (Catalytic Wall BCs, 2026-09)

有限速率催化 (Scott γ-碰撞模型) 与超催化壁面, 以**新增边界条件 + 求解器侧能量完成项**实现,
不修改 `YEqn.H` 的物种方程结构。计算刚性问题暂不处理 (现阶段只求数学物理过程正确与软件适配)。

#### 物理模型 (γ-collision form)

表面复合速率 (E-R 饱和覆盖, 无因子 2, 无 Arrhenius):
- ω = γ·ρ_w·Y_w·√(R_s·T_w/(2π)) = k_c·ρ_w·Y_w,  k_c = γ·√(R_s·T_w/(2π)) [m/s]
- R_s = RR/W, RR = 1000·R [J/(kmol·K)], W 为原子摩尔质量 [kg/kmol] (BC 字典 `molarMass`)

壁面原子质量平衡 −ρD·∂Y/∂n = k_c·ρ·Y_w 以 mixed (Robin) 形式离散:

    valueFraction = k_c·ρ / (k_c·ρ + ρD·δ_c),  refValue = 0,  refGrad = 0

#### 三个新边界条件 (libhyTwoFoam)

| 类型 | 基类 | 用途 |
|---|---|---|
| `nonEqScottCatalyticY` | mixed | 原子组分 (如 N): 有限速率复合 Robin 条件; 关键字 `catalyticEfficiency` (γ), `molarMass`, 可选 `Twall` |
| `nonEqScottCatalyticProductY` | mixed | 复合产物 (如 N2): 梯度条件 snGrad = ω/ρD (`atomField` 指向原子场); 与原子 BC 一起闭 合壁面质量平衡 (ΣIs = 0 严格成立) |
| `nonEqSuperCatalyticY` | fixedValue | 超催化: 直接给定壁面组分 (γ→∞ 供给受限极限, 与 γ=1 的动力学上限相区别); 壁面 不穿透由 Stefan 修正自动保证 |

惰性组分 (如 NO、Ar) 用 `zeroGradient`; 纯氮冻结流可继续用 `nonCatalytic` (等价 zeroGradient)。

#### 壁面温度同源

催化与温度跳跃模型共用同一壁面温度, 按优先级解析:
1. BC 字典 `Twall` (均匀或逐面);
2. Tt 为 `nonEqSmoluchowskiJumpT` 时取其 `Twall()` (需 recompile 后的访问器);
3. Tt 为 `fixedValue` 时取其值;
4. 否则用气侧 Tt 面值并警告一次。

#### 能量完成项 (求解器侧)

`multiSpeciesHeatSource()`/`VEHeatSource()` 只填充内部单元 (边界面值为零), 催化组分通量输运
的壁面焓通量因此从未进入能量方程。新增 `eqns/catalyticWallSources.H` (include 于 `YEqn.H`
末尾) 逐迭代构造完成通量场的边界面值:

- Is = −ρD·snGrad(Y) (与 Y 方程 laplacian 同一 ρD 场), Js = Is − Y_w·ΣIs (Stefan, 镜像
  `getDiffusiveWallHeatFlux()`)
- **气侧用显焓** (hs, Hvels — hyStrath 的 e 为显能形式, 镜像内部源项): 生成能份额 ω·Δh_f 只进
  壁面报告, 由气相化学池 (组分变化) 承担, 与气相反应处理一致
- **壁面报告用绝对焓** (ha): `wallHeatFlux_cat` 含全部复合反应热 ω·Δh_recomb, 指向壁面为正
- 能量模式分配: 单 T → `catalyticFluxE` 进 e 方程; 单 Tv → `catalyticFluxTr`/`FluxVe` 进
  et/ev 方程; 逐分子 Tv → `catalyticFluxTr` + 逐物种 `catalyticFluxVeK`, 与 `eEqnViscous.H`
  的三路分配逐项镜像
- 算符符号与内部项一致 (`eqn += fvc::div(F)`): 气相失去复合原子的显焓, ev 池获得产物分子
  在壁面温度下的振动能

输出: `wallHeatFlux_cat` 场 (NO_WRITE, write.H 显式写) + 日志逐壁面积分催化热流与占总热流
百分比 (Goulard 型诊断)。

#### 配置示例 (air-5, N/O 复合)

```cpp
// 0/N
wall { type nonEqScottCatalyticY; catalyticEfficiency 0.1; molarMass 14.0067; value uniform 0.1; }
// 0/O
wall { type nonEqScottCatalyticY; catalyticEfficiency 0.05; molarMass 15.999; value uniform 0.05; }
// 0/N2
wall { type nonEqScottCatalyticProductY; atomField N; value uniform 0.79; }
// 0/O2
wall { type nonEqScottCatalyticProductY; atomField O; value uniform 0.2; }
// 0/NO, 0/Ar: zeroGradient (惰性)
```

超催化 (完全复合壁) 则全部组分用 `nonEqSuperCatalyticY` 给定平衡组分 (如 N2: 1, 其余: 0)。

#### 回归不变量

- 无催化 BC → `catalyticWallPresent = false` → 求解器侧全部改动跳过 (逐位不变)
- γ = 0 → valueFraction = 0 → 精确 zeroGradient (逐位不变)
- 壁面检测仅按类型名 (纯 word 比较), 求解器无 BC 库链接期依赖
- 启动校验: 催化需多组分扩散模型 + 求解组分方程; 同壁面 Scott 与 Super 互斥 (FatalError);
  Scott 无 Product 配套 → 警告 (质量不闭合)

#### 已知近似与限制 (现阶段)

- Robin 条件只作用于 Fick 份额 −ρD·∂Y/∂n; 压力/热扩散 (JGradp/JGradT) 在壁面的微小份额未
  折入动力学平衡 (完成项与报告同样只取 Fick+Stefan 份额, 内部自洽; SU2 同样如此处理)。
- 产物 BC 为通量注入型梯度条件, 强注入+粗网格时 Y_w,产物 可瞬时略超 1 (自校正, 不影响通量
  精度; 未做裁剪)。
- 上游 `multiSpeciesTransportModel::getDiffusiveWallHeatFlux()` 的热流累加存在缺陷
  (`multiSpeciesTransportModel.C:644` 用 `=` 而非 `+=`, 报告的 diffusive 项只含最后一个
  重组分); `wallHeatFlux_cat` 按全重组分正确求和, 不受影响; 总热流占比因此可能偏大。
  修复该上游缺陷需另行确认。
- inviscid 模式下催化对能量收支无作用 (启动时警告)。
- 刚性 (Da→∞ 时 Robin 退化 Dirichlet 的隐性矩阵性态) 未处理, 由时间步控制承担。
