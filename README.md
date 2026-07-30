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

### 通量格式

在 `system/fvSchemes` 中通过 `fluxScheme` 关键字选择:

| fluxScheme | 精度 | 说明 |
|------------|------|------|
| `Kurganov` | 2 阶 | 原始 Kurganov 中心格式 (默认) |
| `Tadmor` | 2 阶 | 原始 Tadmor 中心格式 |
| `MUSCL_AUSM` | 3 阶 | MUSCL (κ=1/3, van Albada) + AUSM+-up |
| `WENOZ_AUSM` | 5 阶 | WENO-Z (Borges 2008) + AUSM+-up |
| `WENOZ7_AUSM` | 7 阶 | WENO-Z + AUSM+-up |

注: AUSM 系列仅适用于结构化网格 (blockMesh 六面体), 非轴对齐面自动回退为 1 阶 upwind.

### 梯度格式

| gradScheme | 精度 | 说明 |
|------------|------|------|
| `sixth` | 6 阶 | 六阶紧致梯度 (需单独编译到 finiteVolume 库) |

### 新增/修改文件

```
hy2Foam/numerics/
├── AUSM_Flux.H                  # [新] AUSM+-up 通量函数
├── MUSCL_Kernels.H              # [新] MUSCL 插值核
├── WENOZ_Kernels.H              # [新] WENO-Z 插值核
├── highOrderKernels.H           # [新] 模板/重构辅助函数
├── highOrderData.H              # [新] WENO 权重数据
├── MUSCL_WENO_Reconstruct.H     # [新] 每步面值重构
├── structuredStencilBuild.H     # [新] 结构化模板构建
├── readFluxScheme.H             # [改] 注册新通量方案
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

修复后的算子分裂:
- `evEqnViscous`:  rho * d(ev)/dt = div(kappave * grad(Tv))
- `eEqnViscous`:   rho * d(et)/dt = div(kappatr * grad(Ttr))
- 合计:           rho * d(e)/dt  = div(kappa * grad(T))  ← 物理正确

### 原有功能完全保留

当 `fluxScheme` 为 `Kurganov` 或 `Tadmor` 时, 所有新代码无操作, 结果与原始 hyStrath 完全一致.
