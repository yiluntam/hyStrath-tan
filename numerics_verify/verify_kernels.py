# -*- coding: utf-8 -*-
# 数值验证 MUSCL / WENO-Z 核函数 (严格复现 C++ 实现)
#  1) L/R 面状态对光滑解的收敛阶 (期望: MUSCL-3 -> 3, WENOZ5 -> 5, WENOZ7 -> 7)
#  2) beta4JS 与精确符号积分对比 (逐系数核对)
#  3) tau7 组合 (1,3,-3,-1) 的 O(h^7) 论断 vs (1,-3,3,-1) 的 O(h^5)
import math

# ---------- C++ 核函数的 Python 复现 ----------

def musclLimiterPhi(r, ltype):
    if ltype == 0:   # minmod
        return max(min(1.0, r), 0.0)
    if ltype == 1:   # vanLeer
        a = abs(r)
        return (r + a) / (1.0 + a)
    if ltype == 2:   # superbee
        return max(max(min(2.0*r, 1.0), min(r, 2.0)), 0.0)
    # vanAlbada (未钳位, 与 C++ 一致)
    return r*(r + 1.0)/(r*r + 1.0)

def muscl3FaceL(uM1, uC, uP1, kappa, ltype):
    dM = uC - uM1
    dP = uP1 - uC
    phi = musclLimiterPhi(dM/(dP + 1.0e-20), ltype)
    return uC + 0.25*((1.0 - kappa)*dM + (1.0 + kappa)*phi*dP)

def muscl3FaceR(uC, uP1, uP2, kappa, ltype):
    dM = uP1 - uC          # 界面侧
    dP = uP2 - uP1         # 迎风侧
    phi = musclLimiterPhi(dP/(dM + 1.0e-20), ltype)
    return uP1 - 0.25*((1.0 - kappa)*dP + (1.0 + kappa)*phi*dM)

def wenoZ5_faceL(uM2, uM1, uC, uP1, uP2):
    eps = 1e-40
    p0 = (2.0*uM2 - 7.0*uM1 + 11.0*uC)/6.0
    p1 = (-uM1 + 5.0*uC + 2.0*uP1)/6.0
    p2 = (2.0*uC + 5.0*uP1 - uP2)/6.0
    b0 = (13.0/12.0)*(uM2 - 2.0*uM1 + uC)**2 + 0.25*(uM2 - 4.0*uM1 + 3.0*uC)**2
    b1 = (13.0/12.0)*(uM1 - 2.0*uC + uP1)**2 + 0.25*(uM1 - uP1)**2
    b2 = (13.0/12.0)*(uC - 2.0*uP1 + uP2)**2 + 0.25*(3.0*uC - 4.0*uP1 + uP2)**2
    tau5 = abs(b0 - b2)
    a0 = 0.1*(1.0 + (tau5/(b0 + eps))**2)
    a1 = 0.6*(1.0 + (tau5/(b1 + eps))**2)
    a2 = 0.3*(1.0 + (tau5/(b2 + eps))**2)
    s = a0 + a1 + a2
    return (a0*p0 + a1*p1 + a2*p2)/s

def wenoZ5_faceR(uM1, uC, uP1, uP2, uP3):
    eps = 1e-40
    p0 = (2.0*uP3 - 7.0*uP2 + 11.0*uP1)/6.0
    p1 = (-uP2 + 5.0*uP1 + 2.0*uC)/6.0
    p2 = (2.0*uP1 + 5.0*uC - uM1)/6.0
    b0 = (13.0/12.0)*(uP3 - 2.0*uP2 + uP1)**2 + 0.25*(uP3 - 4.0*uP2 + 3.0*uP1)**2
    b1 = (13.0/12.0)*(uP2 - 2.0*uP1 + uC)**2 + 0.25*(uP2 - uC)**2
    b2 = (13.0/12.0)*(uP1 - 2.0*uC + uM1)**2 + 0.25*(3.0*uP1 - 4.0*uC + uM1)**2
    tau5 = abs(b0 - b2)
    a0 = 0.1*(1.0 + (tau5/(b0 + eps))**2)
    a1 = 0.6*(1.0 + (tau5/(b1 + eps))**2)
    a2 = 0.3*(1.0 + (tau5/(b2 + eps))**2)
    s = a0 + a1 + a2
    return (a0*p0 + a1*p1 + a2*p2)/s

def beta4JS(f0, f1, f2, f3, p):
    v1 = f1 - f0
    v2 = f2 - f1
    v3 = f3 - f2
    if p == 0:
        return (2107.0*v1*v1 - 5188.0*v1*v2 + 1854.0*v1*v3
              + 3708.0*v2*v2 - 2788.0*v2*v3 + 547.0*v3*v3)/240.0
    if p == 1:
        return (547.0*v1*v1 - 1428.0*v1*v2 + 494.0*v1*v3
             + 1468.0*v2*v2 - 1108.0*v2*v3 + 267.0*v3*v3)/240.0
    if p == 2:
        return (267.0*v1*v1 - 1108.0*v1*v2 + 494.0*v1*v3
             + 1468.0*v2*v2 - 1428.0*v2*v3 + 547.0*v3*v3)/240.0
    return (547.0*v1*v1 - 2788.0*v1*v2 + 1854.0*v1*v3
         + 3708.0*v2*v2 - 5188.0*v2*v3 + 2107.0*v3*v3)/240.0

def wenoZ7_faceL(uM3, uM2, uM1, uC, uP1, uP2, uP3):
    eps = 1e-40
    v0 = (-3.0*uM3 + 13.0*uM2 - 23.0*uM1 + 25.0*uC)/12.0
    v1 = (uM2 - 5.0*uM1 + 13.0*uC + 3.0*uP1)/12.0
    v2 = (-uM1 + 7.0*uC + 7.0*uP1 - uP2)/12.0
    v3 = (3.0*uC + 13.0*uP1 - 5.0*uP2 + uP3)/12.0
    b0 = beta4JS(uM3, uM2, uM1, uC, 3)
    b1 = beta4JS(uM2, uM1, uC, uP1, 2)
    b2 = beta4JS(uM1, uC, uP1, uP2, 1)
    b3 = beta4JS(uC, uP1, uP2, uP3, 0)
    tau7 = abs(b0 + 3.0*b1 - 3.0*b2 - b3)
    g = (1.0/35.0, 12.0/35.0, 18.0/35.0, 4.0/35.0)
    b = (b0, b1, b2, b3)
    v = (v0, v1, v2, v3)
    a = [g[k]*(1.0 + (tau7/(b[k] + eps))**2) for k in range(4)]
    s = sum(a)
    return sum(a[k]*v[k] for k in range(4))/s

def wenoZ7_faceR(uM2, uM1, uC, uP1, uP2, uP3, uP4):
    eps = 1e-40
    v0 = (-3.0*uP4 + 13.0*uP3 - 23.0*uP2 + 25.0*uP1)/12.0
    v1 = (uP3 - 5.0*uP2 + 13.0*uP1 + 3.0*uC)/12.0
    v2 = (-uP2 + 7.0*uP1 + 7.0*uC - uM1)/12.0
    v3 = (3.0*uP1 + 13.0*uC - 5.0*uM1 + uM2)/12.0
    b0 = beta4JS(uP1, uP2, uP3, uP4, 0)
    b1 = beta4JS(uC, uP1, uP2, uP3, 1)
    b2 = beta4JS(uM1, uC, uP1, uP2, 2)
    b3 = beta4JS(uM2, uM1, uC, uP1, 3)
    tau7 = abs(b0 + 3.0*b1 - 3.0*b2 - b3)
    g = (1.0/35.0, 12.0/35.0, 18.0/35.0, 4.0/35.0)
    b = (b0, b1, b2, b3)
    v = (v0, v1, v2, v3)
    a = [g[k]*(1.0 + (tau7/(b[k] + eps))**2) for k in range(4)]
    s = sum(a)
    return sum(a[k]*v[k] for k in range(4))/s

# ---------- 光滑试验函数与胞平均 ----------

def f_exact(x):
    return math.sin(2.0*math.pi*x) + 0.3*math.cos(6.0*math.pi*x)

def cell_avg(i, N):
    # 10 点 Gauss-Legendre 胞平均, h = 1/N
    h = 1.0/N
    # Gauss-Legendre 节点/权 (n=10) 于 [-1,1]
    t = (-0.973906528517171720, -0.865063366688984511, -0.679409568299024406,
         -0.433395394129247191, -0.148874338981631211, 0.148874338981631211,
         0.433395394129247191, 0.679409568299024406, 0.865063366688984511,
         0.973906528517171720)
    w = (0.066671344308688138, 0.149451349150580593, 0.219086362515982044,
         0.269266719309996355, 0.295524224714752868, 0.295524224714752868,
         0.269266719309996355, 0.219086362515982044, 0.149451349150580593,
         0.066671344308688138)
    a, b = i*h, (i+1)*h
    s = 0.0
    for k in range(10):
        x = 0.5*(b - a)*t[k] + 0.5*(a + b)
        s += w[k]*f_exact(x)
    return 0.5*s*(b - a)/h*h  # = 0.5*(b-a)*sum(w f)/h * h -> 平均

def cell_avg2(i, N):
    h = 1.0/N
    t = (-0.973906528517171720, -0.865063366688984511, -0.679409568299024406,
         -0.433395394129247191, -0.148874338981631211, 0.148874338981631211,
         0.433395394129247191, 0.679409568299024406, 0.865063366688984511,
         0.973906528517171720)
    w = (0.066671344308688138, 0.149451349150580593, 0.219086362515982044,
         0.269266719309996355, 0.295524224714752868, 0.295524224714752868,
         0.269266719309996355, 0.219086362515982044, 0.149451349150580593,
         0.066671344308688138)
    a, b = i*h, (i+1)*h
    s = sum(w[k]*f_exact(0.5*(b-a)*t[k] + 0.5*(a+b)) for k in range(10))
    return s/2.0  # ∫f dx / h = 0.5*sum (因 (b-a)=h)

def order_table(name, err):
    print("  %s" % name)
    prev = None
    for h, (el, er) in err.items():
        if prev is None:
            print("    h=%.4f  errL=%.3e  errR=%.3e" % (h, el, er))
        else:
            pl = math.log(prev[0]/el, 2.0)
            pr = math.log(prev[1]/er, 2.0)
            print("    h=%.4f  errL=%.3e (%.2f)  errR=%.3e (%.2f)"
                  % (h, el, pl, er, pr))
        prev = (el, er)

def run_convergence():
    print("=== 1) 收敛阶验证 (光滑周期解, 内部面) ===")
    results = {}
    for scheme in ("MUSCL3", "WENOZ5", "WENOZ7"):
        errs = {}
        for N in (40, 80, 160, 320):
            h = 1.0/N
            u = [cell_avg2(i, N) for i in range(N+8)]  # 周期延拓索引偏移 +4
            def U(i):  # 周期
                return u[i % (N+8) - 0] if False else u[i]
            # 用周期数据填充: 直接生成 N 个胞并循环索引
            uc = [cell_avg2(i, N) for i in range(N)]
            def UU(i):
                return uc[i % N]
            el = er = 0.0
            nface = 0
            for j in range(4, N-4):  # 内部面 j+1/2, 留足模板
                # C++ 调用约定: 面 O|N, O=j, N=j+1
                if scheme == "MUSCL3":
                    k = 1.0/3.0
                    L = muscl3FaceL(UU(j-1), UU(j), UU(j+1), k, 3)
                    R = muscl3FaceR(UU(j), UU(j+1), UU(j+2), k, 3)
                elif scheme == "WENOZ5":
                    L = wenoZ5_faceL(UU(j-2), UU(j-1), UU(j), UU(j+1), UU(j+2))
                    R = wenoZ5_faceR(UU(j-1), UU(j), UU(j+1), UU(j+2), UU(j+3))
                else:
                    L = wenoZ7_faceL(UU(j-3), UU(j-2), UU(j-1), UU(j),
                                     UU(j+1), UU(j+2), UU(j+3))
                    R = wenoZ7_faceR(UU(j-2), UU(j-1), UU(j), UU(j+1),
                                     UU(j+2), UU(j+3), UU(j+4))
                exact = f_exact((j+1.0)*h)  # 胞 i 跨 [ih,(i+1)h], 面 j|j+1 在 (j+1)h
                el += (L - exact)**2
                er += (R - exact)**2
                nface += 1
            errs[h] = (math.sqrt(el/nface), math.sqrt(er/nface))
        results[scheme] = errs
        order_table(scheme, errs)
    return results

# ---------- 2) beta4JS vs 精确符号积分 ----------

def beta_exact(f0, f1, f2, f3, p):
    # JS: beta = sum_{l=1..3} h^{2l-1} \int_{cell p} (q^{(l)})^2 dx, h=1
    # q = 过 4 个胞平均的三次多项式 (均匀网格, 胞中心 0,1,2,3)
    # 用幂基在胞中心处展开: q(x) = c0 + c1 x + c2 x^2 + c3 x^3
    # 胞平均 of x^m over cell [k-1/2, k+1/2]:
    #   m=0: 1
    #   m=1: k
    #   m=2: k^2 + 1/12
    #   m=3: k^3 + k/4
    ks = (0, 1, 2, 3)
    fs = (f0, f1, f2, f3)
    # 组装 Vandermode-like 方程 (胞平均), 解 c
    import itertools
    A = []
    for k in ks:
        A.append([1.0, k, k*k + 1.0/12.0, k*k*k + k/4.0])
    # 高斯消元 4x4
    M = [row[:] + [fs[i]] for i, row in enumerate(A)]
    for col in range(4):
        piv = max(range(col, 4), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        for r in range(4):
            if r != col and abs(M[r][col]) > 1e-300:
                fac = M[r][col]/M[col][col]
                for cc in range(col, 5):
                    M[r][cc] -= fac*M[col][cc]
    c = [M[i][4]/M[i][i] for i in range(4)]
    c1, c2, c3 = c[1], c[2], c[3]
    # q' = c1 + 2c2 x + 3c3 x^2 ; q'' = 2c2 + 6c3 x ; q''' = 6c3
    # 积分区: 胞 p: [p-1/2, p+1/2]
    a, b = p - 0.5, p + 0.5
    def int2(poly):  # poly 系数 (最低次在前), 积分平方
        n = len(poly)
        coef = [0.0]*(2*n - 1)
        for i in range(n):
            for j in range(n):
                coef[i + j] += poly[i]*poly[j]
        # 积分 a..b
        s = 0.0
        for i, ci in enumerate(coef):
            s += ci*(b**(i+1) - a**(i+1))/(i+1)
        return s
    qp = [c1, 2*c2, 3*c3]
    qpp = [2*c2, 6*c3]
    qppp = [6*c3]
    return int2(qp) + int2(qpp) + int2(qppp)

def run_beta_check():
    print("=== 2) beta4JS vs 精确符号积分 ===")
    maxdiff = 0.0
    cases = [
        (1.0, -2.0, 3.5, 0.2), (5.0, 5.0, 5.0, 5.0),
        (0.1, 100.0, -3.0, 42.0), (1.0, 2.0, 3.0, 4.0),
        (-7.7, 0.3, 12.0, -0.9),
    ]
    for cs in cases:
        for p in range(4):
            a = beta4JS(*cs, p)
            b = beta_exact(*cs, p)
            d = abs(a - b)/max(abs(b), 1e-12)
            maxdiff = max(maxdiff, d)
            if d > 1e-10:
                print("  MISMATCH cs=%s p=%d  code=%.12g exact=%.12g rel=%.2e"
                      % (cs, p, a, b, d))
    print("  最大相对偏差: %.3e  %s" % (maxdiff,
          "(全部一致)" if maxdiff < 1e-10 else "(存在偏差!)"))

# ---------- 3) tau7 组合量阶 ----------

def run_tau7_scaling():
    print("=== 3) tau7 组合量阶 (光滑解, 系数 O(h^p)) ===")
    for combo, name in (((1, 3, -3, -1), "代码 (1,3,-3,-1)"),
                        ((1, -3, 3, -1), "三阶差分 (1,-3,3,-1)")):
        print("  %s:" % name)
        prev = None
        for N in (20, 40, 80, 160):
            h = 1.0/N
            uc = [cell_avg2(i, N) for i in range(N)]
            def UU(i):
                return uc[i % N]
            s = 0.0
            cnt = 0
            for j in range(4, N-4):
                b0 = beta4JS(UU(j-3), UU(j-2), UU(j-1), UU(j), 3)
                b1 = beta4JS(UU(j-2), UU(j-1), UU(j), UU(j+1), 2)
                b2 = beta4JS(UU(j-1), UU(j), UU(j+1), UU(j+2), 1)
                b3 = beta4JS(UU(j), UU(j+1), UU(j+2), UU(j+3), 0)
                val = abs(combo[0]*b0 + combo[1]*b1 + combo[2]*b2 + combo[3]*b3)
                s += val
                cnt += 1
            mean = s/cnt
            if prev is not None:
                print("    h=%.5f  <tau7>=%.3e  阶=%.2f"
                      % (h, mean, math.log(prev/mean, 2.0)))
            else:
                print("    h=%.5f  <tau7>=%.3e" % (h, mean))
            prev = mean

# ---------- 4) 镜像对称性 ----------

def run_mirror():
    print("=== 4) 镜像对称性: R(反射数据) == L(原数据) ===")
    import random
    random.seed(42)
    N = 200
    uc = [cell_avg2(i, N) for i in range(N)]
    def UU(i):
        return uc[i % N]
    d5 = d7 = 0.0
    for j in range(4, N-4):
        # 恒等式: R(d) == L(d∘refl), refl: k -> 2j+1-k (关于面 j|j+1 反射)
        # R 状态参数表为 (j-1..j+3), 反射后 = L 的参数表 (j+3..j-1)
        L = wenoZ5_faceL(UU(j+3), UU(j+2), UU(j+1), UU(j), UU(j-1))
        Rr = wenoZ5_faceR(UU(j-1), UU(j), UU(j+1), UU(j+2), UU(j+3))
        d5 = max(d5, abs(L - Rr))
        L7 = wenoZ7_faceL(UU(j+4), UU(j+3), UU(j+2), UU(j+1),
                          UU(j), UU(j-1), UU(j-2))
        Rr7 = wenoZ7_faceR(UU(j-2), UU(j-1), UU(j), UU(j+1),
                           UU(j+2), UU(j+3), UU(j+4))
        d7 = max(d7, abs(L7 - Rr7))
    print("  WENOZ5 max|L - R(refl)| = %.3e" % d5)
    print("  WENOZ7 max|L - R(refl)| = %.3e" % d7)

if __name__ == "__main__":
    run_beta_check()
    run_tau7_scaling()
    run_mirror()
    run_convergence()
