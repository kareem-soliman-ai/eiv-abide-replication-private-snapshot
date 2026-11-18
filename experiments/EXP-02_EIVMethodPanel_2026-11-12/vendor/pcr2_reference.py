"""
PCR-2 reference derivation check: clouds (per-point measurement covariance)
around a hazy line (intrinsic scatter tau^2).

Everything here is deliberately small and self-contained so each equation in
the write-up can be checked against a number.

Model, per point i:
    xi_i   ~ G                      latent true x   (structural: N(mu, T))
    eta_i  = a + b*xi_i + eps_i,    eps_i ~ N(0, tau^2)        ("haze")
    (x_i, y_i) = (xi_i, eta_i) + w_i,  w_i ~ N(0, Sigma_i)     ("cloud")
    Sigma_i = [[sx_i^2, rho_i sx_i sy_i], [rho_i sx_i sy_i, sy_i^2]]
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, minimize

rng_global = np.random.default_rng(20260927)


# ----------------------------------------------------------------------------
# Estimators
# ----------------------------------------------------------------------------
def ols(x, y):
    b = np.cov(x, y, bias=True)[0, 1] / np.var(x)
    return b, y.mean() - b * x.mean()


def regression_calibration(x, y, sx):
    """Moment / disattenuation estimator: b = Sxy / (Sxx - mean sx^2)."""
    sxx = np.var(x) - np.mean(sx**2)
    b = np.cov(x, y, bias=True)[0, 1] / sxx
    return b, y.mean() - b * x.mean()


def york(x, y, sx, sy, tau2=0.0, rho=None, b0=None, tol=1e-12, maxit=500):
    """York et al. (2004) iteration, with the y-error variance augmented by tau2.

    Minimises  sum_i r_i^2 / V_i,   V_i = tau2 + sy_i^2 + b^2 sx_i^2 - 2 b rho_i sx_i sy_i
    """
    n = len(x)
    rho = np.zeros(n) if rho is None else np.asarray(rho, float)
    sy_tot = np.sqrt(sy**2 + tau2)
    r_eff = rho * sy / sy_tot            # correlation between x-error and total y-deviation
    wx, wy = 1.0 / sx**2, 1.0 / sy_tot**2
    alpha = np.sqrt(wx * wy)
    b = ols(x, y)[0] if b0 is None else b0
    for it in range(1, maxit + 1):
        W = wx * wy / (wx + b * b * wy - 2 * b * r_eff * alpha)
        Xbar, Ybar = np.sum(W * x) / W.sum(), np.sum(W * y) / W.sum()
        U, V = x - Xbar, y - Ybar
        beta = W * (U / wy + b * V / wx - (b * U + V) * r_eff / alpha)
        b_new = np.sum(W * beta * V) / np.sum(W * beta * U)
        if abs(b_new - b) < tol * max(1.0, abs(b)):
            b = b_new
            break
        b = b_new
    W = wx * wy / (wx + b * b * wy - 2 * b * r_eff * alpha)
    a = np.sum(W * y) / W.sum() - b * np.sum(W * x) / W.sum()
    chi2 = float(np.sum(W * (y - a - b * x) ** 2))
    return b, a, chi2, it


def hazy_york(x, y, sx, sy, rho=None):
    """York with tau2 chosen so that the reduced chi-square equals one.

    This is the astrophysics intrinsic-scatter convention (Tremaine et al. 2002)
    and a method-of-moments estimator of tau2: E[r_i^2] = V_i at the truth.
    """
    n = len(x)

    def g(t2):
        return york(x, y, sx, sy, tau2=t2, rho=rho)[2] - (n - 2)

    if g(0.0) <= 0:
        t2 = 0.0
    else:
        hi = 10.0 * np.var(y)
        while g(hi) > 0:
            hi *= 10
        t2 = brentq(g, 0.0, hi, xtol=1e-10)
    b, a, chi2, it = york(x, y, sx, sy, tau2=t2, rho=rho)
    return b, a, t2


def structural_ml(x, y, sx, sy, rho=None):
    """Gaussian structural EIV with haze: z_i ~ N(m, Lambda + Sigma_i).

    theta = (mu, log T, a, b, tau2).  Direct maximisation, no alternating step.
    """
    n = len(x)
    rho = np.zeros(n) if rho is None else np.asarray(rho, float)
    sxx, syy, sxy = sx**2, sy**2, rho * sx * sy

    def nll(theta):
        mu, logT, a, b, tau2 = theta
        T = np.exp(logT)
        c11 = T + sxx
        c12 = b * T + sxy
        c22 = b * b * T + tau2 + syy
        det = c11 * c22 - c12**2
        dx, dy = x - mu, y - (a + b * mu)
        quad = (c22 * dx * dx - 2 * c12 * dx * dy + c11 * dy * dy) / det
        return 0.5 * np.mean(np.log(det) + quad)

    T0 = max(np.var(x) - np.mean(sxx), 0.05)
    b0 = np.cov(x, y, bias=True)[0, 1] / T0
    a0 = y.mean() - b0 * x.mean()
    tau0 = max(np.var(y) - b0**2 * T0 - np.mean(syy), 0.05)
    best = None
    for bstart in (b0, 0.0):
        res = minimize(nll, [x.mean(), np.log(T0), a0, bstart, tau0],
                       method="L-BFGS-B",
                       bounds=[(None, None), (np.log(1e-8), np.log(1e4)),
                               (None, None), (-50, 50), (0.0, 1e4)])
        if best is None or res.fun < best.fun:
            best = res
    mu, logT, a, b, tau2 = best.x
    return b, a, tau2, np.exp(logT), mu


def old_pcr_rule(x, y, sx, sy, tol=1e-8, maxit=500):
    """The original PCR alternating rule (pcr.py): projection then WLS with
    combined-variance weights 1/(sy^2 + b^2 sx^2)."""
    b, a = ols(x, y)
    for it in range(1, maxit + 1):
        den = sy**2 + b**2 * sx**2
        xh = (sy**2 * x + sx**2 * b * (y - a)) / den
        w = 1.0 / den
        xb, yb = np.sum(w * xh) / w.sum(), np.sum(w * y) / w.sum()
        b_new = np.sum(w * (xh - xb) * (y - yb)) / np.sum(w * (xh - xb) ** 2)
        a_new = yb - b_new * xb
        change = max(abs(b_new - b), abs(a_new - a))
        b, a = b_new, a_new
        if not np.isfinite(change):
            return b, a, it, False
        if change < tol:
            return b, a, it, True
    return b, a, it, False


def coordinate_descent(x, y, sx, sy, tol=1e-8, maxit=5000):
    """Objective-coherent alternating minimisation of the joint Mahalanobis
    objective: projection step, then WLS of y on x_hat with weights 1/sy^2."""
    b, a = ols(x, y)
    for it in range(1, maxit + 1):
        xh = x + b * sx**2 * (y - a - b * x) / (sy**2 + b**2 * sx**2)
        w = 1.0 / sy**2
        xb, yb = np.sum(w * xh) / w.sum(), np.sum(w * y) / w.sum()
        b_new = np.sum(w * (xh - xb) * (y - yb)) / np.sum(w * (xh - xb) ** 2)
        a_new = yb - b_new * xb
        change = max(abs(b_new - b), abs(a_new - a))
        b, a = b_new, a_new
        if change < tol:
            return b, a, it, True
    return b, a, it, False


# ----------------------------------------------------------------------------
# Per-point adjustment (the "how much per point" answer)
# ----------------------------------------------------------------------------
def adjust_functional_haze(x, y, sx, sy, a, b, tau2, rho=None):
    """u_hat_i = z_i - Sigma_i c * r_i / V_i, with c = (-b, 1),
    V_i = tau2 + c' Sigma_i c.  Conditions on the residual only (flat prior on xi).
    Returns adjusted (xi_hat, eta_hat), measurement share f_i, and residual retained."""
    n = len(x)
    rho = np.zeros(n) if rho is None else np.asarray(rho, float)
    sxy = rho * sx * sy
    r = y - a - b * x
    cSc = sy**2 + b * b * sx**2 - 2 * b * sxy          # c' Sigma c
    V = tau2 + cSc
    # Sigma c = (sx^2*(-b) + sxy*1, sxy*(-b) + sy^2)
    Sc_x = -b * sx**2 + sxy
    Sc_y = -b * sxy + sy**2
    xi_hat = x - Sc_x * r / V
    eta_hat = y - Sc_y * r / V
    f = cSc / V                                          # fraction of residual removed
    return xi_hat, eta_hat, f, r * tau2 / V


def adjust_structural(x, y, sx, sy, a, b, tau2, T, mu, rho=None):
    """Full Kalman-gain posterior mean under the structural prior xi ~ N(mu, T)."""
    n = len(x)
    rho = np.zeros(n) if rho is None else np.asarray(rho, float)
    m = np.array([mu, a + b * mu])
    Lam = np.array([[T, b * T], [b * T, b * b * T + tau2]])
    out = np.empty((n, 2))
    for i in range(n):
        Sig = np.array([[sx[i]**2, rho[i] * sx[i] * sy[i]],
                        [rho[i] * sx[i] * sy[i], sy[i]**2]])
        K = Lam @ np.linalg.inv(Lam + Sig)
        out[i] = m + K @ (np.array([x[i], y[i]]) - m)
    return out[:, 0], out[:, 1]


# ----------------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------------
def simulate(n, tau, a=1.0, b=2.0, rng=rng_global, scales=((0.3, 0.5), (0.9, 1.0)),
             xi_sd=1.0):
    xi = rng.normal(0, xi_sd, n)
    grp = rng.integers(0, len(scales), n)
    sx = np.array([scales[g][0] for g in grp])
    sy = np.array([scales[g][1] for g in grp])
    eta = a + b * xi + rng.normal(0, tau, n)
    x = xi + rng.normal(0, sx)
    y = eta + rng.normal(0, sy)
    return x, y, sx, sy, xi, eta


def rmse(v, truth):
    v = np.asarray(v, float)
    return float(np.sqrt(np.mean((v - truth) ** 2)))


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    B_TRUE = 2.0
    print("=" * 78)
    print("1. Intrinsic-scatter cells (n=500, 100 reps, sigma groups (.3,.5)/(.9,1.0))")
    print("=" * 78)
    print(f"{'tau':>4} {'OLS':>8} {'RC':>8} {'York t2=0':>10} {'oldPCR':>8} {'(conv)':>7} "
          f"{'hazyYork':>9} {'structML':>9} {'tau2_hat(hY)':>12} {'tau2_hat(ML)':>12}")
    for tau in (0.0, 1.0, 3.0):
        res = {k: [] for k in ("ols", "rc", "york", "pcr", "hy", "ml")}
        conv = 0
        t2s_hy, t2s_ml = [], []
        for rep in range(100):
            x, y, sx, sy, xi, eta = simulate(500, tau)
            res["ols"].append(ols(x, y)[0])
            res["rc"].append(regression_calibration(x, y, sx)[0])
            res["york"].append(york(x, y, sx, sy)[0])
            bp, ap, itp, cp = old_pcr_rule(x, y, sx, sy)
            if cp:
                res["pcr"].append(bp)
                conv += 1
            bh, ah, t2h = hazy_york(x, y, sx, sy)
            res["hy"].append(bh)
            t2s_hy.append(t2h)
            bm, am, t2m, Tm, mum = structural_ml(x, y, sx, sy)
            res["ml"].append(bm)
            t2s_ml.append(t2m)
        print(f"{tau:4.1f} {rmse(res['ols'], B_TRUE):8.4f} {rmse(res['rc'], B_TRUE):8.4f} "
              f"{rmse(res['york'], B_TRUE):10.4f} "
              f"{(rmse(res['pcr'], B_TRUE) if res['pcr'] else float('nan')):8.4f} {conv:6d}% "
              f"{rmse(res['hy'], B_TRUE):9.4f} {rmse(res['ml'], B_TRUE):9.4f} "
              f"{np.mean(t2s_hy):12.4f} {np.mean(t2s_ml):12.4f}   (true tau2={tau**2:.1f})")

    print()
    print("=" * 78)
    print("2. NHANES-like stress cell: correlation ~0.02 (n=500, tau=1, b=0.03, sx=.3, sy=.2)")
    print("=" * 78)
    res = {k: [] for k in ("ols", "york", "hy", "ml", "rc")}
    chi2nu = []
    for rep in range(100):
        x, y, sx, sy, xi, eta = simulate(500, 1.0, a=0.0, b=0.03,
                                         scales=((0.3, 0.2),))
        res["ols"].append(ols(x, y)[0])
        by, ay, chi2, _ = york(x, y, sx, sy)
        res["york"].append(by)
        chi2nu.append(chi2 / (len(x) - 2))
        res["hy"].append(hazy_york(x, y, sx, sy)[0])
        res["ml"].append(structural_ml(x, y, sx, sy)[0])
        res["rc"].append(regression_calibration(x, y, sx)[0])
    print(f"true b = 0.03; mean reduced chi2 under tau2=0: {np.mean(chi2nu):.2f}")
    for k in res:
        v = np.array(res[k])
        print(f"  {k:5s} mean={v.mean():8.4f}  RMSE={rmse(v, 0.03):8.4f}")

    print()
    print("=" * 78)
    print("3. Per-point adjustment: fractions differ per point; tau2=0 corner == old A-step")
    print("=" * 78)
    x, y, sx, sy, xi, eta = simulate(500, 1.0)
    bm, am, t2m, Tm, mum = structural_ml(x, y, sx, sy)
    print(f"structural ML: b={bm:.4f} a={am:.4f} tau2={t2m:.4f} (true 1.0) T={Tm:.4f} mu={mum:.4f}")
    xh, yh, f, kept = adjust_functional_haze(x, y, sx, sy, am, bm, t2m)
    xs, ys = adjust_structural(x, y, sx, sy, am, bm, t2m, Tm, mum)
    print(f"{'i':>3} {'sx':>4} {'sy':>4} {'resid r':>8} {'f_i removed':>11} {'dx':>7} {'dy':>7} "
          f"{'resid after':>11} | {'dx(struct)':>10} {'dy(struct)':>10}")
    for i in range(8):
        print(f"{i:3d} {sx[i]:4.1f} {sy[i]:4.1f} {y[i]-am-bm*x[i]:8.3f} {f[i]:11.3f} "
              f"{xh[i]-x[i]:7.3f} {yh[i]-y[i]:7.3f} {yh[i]-am-bm*xh[i]:11.3f} | "
              f"{xs[i]-x[i]:10.3f} {ys[i]-y[i]:10.3f}")
    # Does adjustment improve recovery of the latent location?
    print(f"RMSE to latent (xi,eta): raw obs  x:{rmse(x, xi):.3f}  y:{rmse(y, eta):.3f}")
    print(f"                         functional-haze  x:{rmse(xh, xi):.3f}  y:{rmse(yh, eta):.3f}")
    print(f"                         structural(Kalman) x:{rmse(xs, xi):.3f}  y:{rmse(ys, eta):.3f}")
    # Full projection (old A-step, tau2=0) for comparison
    xh0, yh0, f0, _ = adjust_functional_haze(x, y, sx, sy, am, bm, 0.0)
    print(f"                         full projection (tau2=0) x:{rmse(xh0, xi):.3f}  y:{rmse(yh0, eta):.3f}")
    # Corner check: tau2=0 adjustment equals the old A-step exactly
    den = sy**2 + bm**2 * sx**2
    xh_old = (sy**2 * x + sx**2 * bm * (y - am)) / den
    print(f"max |x_hat(tau2=0) - old A-step x_hat| = {np.max(np.abs(xh0 - xh_old)):.2e}")

    print()
    print("=" * 78)
    print("4. Why the alternating rule stalls: iterations vs kappa = b^2 sx^2/(sy^2 + b^2 sx^2)")
    print("=" * 78)
    print(f"{'sx':>4} {'kappa':>6} {'coordDesc iters':>15} {'conv':>5} {'oldPCR iters':>13} {'conv':>5} {'York iters':>10}")
    for sxv in (0.1, 0.5, 1.0, 2.0, 3.0):
        x, y, sx, sy, xi, eta = simulate(500, 0.0, scales=((sxv, 0.5),))
        kappa = B_TRUE**2 * sxv**2 / (0.25 + B_TRUE**2 * sxv**2)
        _, _, itc, cc = coordinate_descent(x, y, sx, sy)
        _, _, itp, cp = old_pcr_rule(x, y, sx, sy)
        _, _, _, ity = york(x, y, sx, sy)
        print(f"{sxv:4.1f} {kappa:6.3f} {itc:15d} {str(cc):>5} {itp:13d} {str(cp):>5} {ity:10d}")
