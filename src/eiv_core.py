"""Shared, scale-normalised Gaussian errors-in-variables fits.

Functional fits minimise the Mahalanobis profile sum (no log projected
variance). Structural fits maximise a proper joint Gaussian density for
observations: latent covariance plus known measurement covariance. The
structural latent covariance includes intrinsic conditional outcome scatter.
Rows are assumed independent; measurement scales are supplied, not inferred.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize

METHOD_VERSION = "rasad-eiv/2"


@dataclass
class EIVFit:
    coef: np.ndarray
    intercept: float
    latent: np.ndarray
    posterior_covariances: np.ndarray
    intrinsic_variance: float
    latent_covariance: np.ndarray | None
    latent_mean: np.ndarray
    converged: bool
    iterations: int
    objective_history: list[float]
    log_likelihood: float
    diagnostics: dict


def measurement_covariances(X, y, sigma_X=None, sigma_y=None, covariances=None):
    """Validate inputs; support exact X and positive semidefinite covariance."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    if X.ndim == 1:
        X = X[:, None]
    if X.ndim != 2 or y.ndim != 1 or len(y) != len(X) or len(y) < 3:
        raise ValueError("Require X (n,p), y (n,) and at least three observations")
    n, p = X.shape
    if covariances is None:
        if sigma_X is None or sigma_y is None:
            raise ValueError("Provide sigma_X/sigma_y or observation_covariances")
        sx = np.asarray(sigma_X, float)
        if sx.ndim == 1 and sx.size == n and sx.size != p:
            sx = sx[:, None]
        sx = np.broadcast_to(sx, X.shape)
        sy = np.broadcast_to(np.asarray(sigma_y, float), y.shape)
        if np.any(sx < 0) or np.any(sy <= 0):
            raise ValueError("Require sigma_X >= 0 and sigma_y > 0")
        covariances = np.zeros((n, p+1, p+1))
        diag = np.arange(p+1)
        covariances[:, diag, diag] = np.column_stack([sx*sx, sy*sy])
    C = np.array(covariances, dtype=float, copy=True)
    if C.shape != (n, p+1, p+1):
        raise ValueError(f"observation_covariances must have shape {(n, p+1, p+1)}")
    if not all(np.all(np.isfinite(a)) for a in (X, y, C)):
        raise ValueError("Measurements and covariance must be finite")
    # Test PSD after scaling, so the test does not depend on measurement units.
    scales = np.sqrt(np.maximum(np.diagonal(C, axis1=1, axis2=2), np.finfo(float).tiny))
    correlations = C / scales[:, :, None] / scales[:, None, :]
    if (not np.allclose(correlations, correlations.transpose(0, 2, 1), atol=1e-9)
            or np.any(np.linalg.eigvalsh(correlations) < -1e-9)
            or np.any(C[:, -1, -1] <= 0)):
        raise ValueError("Require symmetric positive semidefinite measurement covariance and positive outcome variance")
    return X, y, C


def fit_eiv(X, y, covariance, *, method="functional", max_iter=200, tol=1e-8, fixed_intrinsic_variance=None):
    if method not in ("functional", "structural"):
        raise ValueError("method must be functional, structural or explicit legacy")
    if max_iter < 1 or tol <= 0:
        raise ValueError("Require max_iter >= 1 and tol > 0")
    z = np.column_stack([X, y])
    n, d = z.shape
    p = d-1
    centre = z.mean(0)
    scale = z.std(0)
    if np.any(scale <= np.finfo(float).eps):
        raise ValueError("Require varying predictors and outcome")
    zn = (z-centre)/scale
    C = covariance/scale[:, None]/scale[None, :]
    # Broadcasting above: scale[:,None] scales rows of each covariance.
    if method == "functional":
        result = _functional(zn, C, max_iter, tol)
    else:
        result = _structural(zn, C, max_iter, tol, None if fixed_intrinsic_variance is None else fixed_intrinsic_variance/scale[-1]**2)
    coef_n, a_n, latent_n, P_n, tau_n, S_n, mu_n, opt, history, diagnostics = result
    coef = coef_n*scale[-1]/scale[:p]
    intercept = float(centre[-1]+scale[-1]*a_n-coef@centre[:p])
    latent = centre + latent_n*scale
    posterior = P_n*scale[:, None]*scale[None, :]
    S = None if S_n is None else S_n*scale[:, None]*scale[None, :]
    mean = centre+mu_n*scale
    log_likelihood = -float(opt.fun)*n
    if method == "structural":
        log_likelihood -= n*(np.log(scale).sum()+d*.5*np.log(2*np.pi))
    diagnostics.update(method=method, method_version=METHOD_VERSION,
                       status="converged" if diagnostics["converged"] else "not_converged",
                       n_observations=n, n_predictors=p,
                       assumptions="independent rows; fixed calibrated measurement covariance",
                       message=str(opt.message))
    return EIVFit(coef, intercept, latent, posterior, float(tau_n*scale[-1]**2),
                  S, mean, diagnostics["converged"], int(opt.nit),
                  [float(v) for v in history], log_likelihood, diagnostics)


def _functional(z, C, max_iter, tol):
    n, d = z.shape
    X, y = z[:, :-1], z[:, -1]
    def objective(theta):
        a, b = theta[0], theta[1:]
        c = np.r_[-b, 1.]
        residual = y-a-X@b
        Sc = C@c
        variance = np.einsum("j,nj->n", c, Sc)
        if np.any(variance <= 0):
            return np.inf, np.zeros(d)
        value = .5*np.mean(residual**2/variance)
        gradient = np.r_[-np.mean(residual/variance),
            np.mean(-X*residual[:, None]/variance[:, None]
                    + Sc[:, :-1]*(residual**2/variance**2)[:, None], axis=0)]
        return float(value), gradient
    start = np.linalg.lstsq(np.column_stack([np.ones(n), X]), y, rcond=None)[0]
    starts = [start]
    moment = z.T@z/n-C.mean(0)
    if np.linalg.eigvalsh(moment[:-1, :-1]).min() > 1e-8:
        starts.append(np.r_[0., np.linalg.solve(moment[:-1, :-1], moment[:-1, -1])])
    fits = []
    for initial in starts:
        history = [objective(initial)[0]]
        opt = minimize(objective, initial, jac=True, method="L-BFGS-B",
                       callback=lambda x: history.append(objective(x)[0]),
                       options={"maxiter":max_iter, "ftol":min(tol*.01, 1e-12),
                                "gtol":min(tol, 1e-7), "maxls":40})
        fits.append((opt, history))
    opt, history = min(fits, key=lambda f:f[0].fun)
    a, b = opt.x[0], opt.x[1:]
    c = np.r_[-b, 1.]
    Sc = C@c
    variance = np.einsum("j,nj->n", c, Sc)
    latent = z-Sc*((y-a-X@b)/variance)[:, None]
    posterior = C-Sc[:, :, None]*Sc[:, None, :]/variance[:, None, None]
    converged = bool(opt.success and np.isfinite(opt.fun) and np.max(np.abs(opt.jac)) < 1e-4)
    diagnostics = dict(converged=converged, objective="functional_profile_chi_square",
                       projected_gradient=float(np.max(np.abs(opt.jac))),
                       starts_succeeded=sum(bool(f[0].success) for f in fits),
                       likelihood_kind="profile objective score; not a marginal Gaussian density")
    return b, a, latent, posterior, 0., None, latent.mean(0), opt, history, diagnostics


def _structural(z, C, max_iter, tol, fixed_tau2=None):
    n, d = z.shape
    rows, cols = np.tril_indices(d)
    diagonal = rows == cols
    moment = np.cov(z, rowvar=False, bias=True)-C.mean(0)
    vals, vecs = np.linalg.eigh(moment)
    initial_cov = (vecs*np.maximum(vals, 1e-5))@vecs.T
    L0 = np.linalg.cholesky(initial_cov)
    if fixed_tau2 is not None:
        if fixed_tau2 < 0:
            raise ValueError("fixed intrinsic variance must be non-negative")
        L0[-1, -1] = np.sqrt(fixed_tau2)
    def unpack(theta):
        L = np.zeros((d, d))
        L[rows, cols] = theta[d:]
        return theta[:d], L
    def objective(theta):
        mu, L = unpack(theta)
        S = L@L.T
        V = S+C
        sign, logdet = np.linalg.slogdet(V)
        if np.any(sign <= 0):
            return np.inf, np.zeros_like(theta)
        inv = np.linalg.inv(V)
        delta = z-mu
        v = np.einsum("nij,nj->ni", inv, delta)
        value = .5*np.mean(logdet+np.einsum("ni,ni->n", delta, v))
        G = .5*np.mean(inv-v[:, :, None]*v[:, None, :], axis=0)
        gradient = np.r_[-v.mean(0), (2*G@L)[rows, cols]]
        return float(value), gradient
    bounds = [(None, None)]*d + [
        ((0. if r == d-1 else 1e-7), None) if r == c else (None, None)
        for r, c in zip(rows, cols)]
    if fixed_tau2 is not None:
        bounds[-1] = (np.sqrt(fixed_tau2), np.sqrt(fixed_tau2))
    starts = [np.r_[z.mean(0), L0[rows, cols]]]
    if fixed_tau2 is None and (np.linalg.eigvalsh(moment).min() <= 0 or not np.allclose(C, C[0])):
        starts.append(np.r_[z.mean(0), np.eye(d)[rows, cols]])
    fits = []
    for initial in starts:
        history = [objective(initial)[0]]
        opt = minimize(objective, initial, jac=True, method="L-BFGS-B", bounds=bounds,
                       callback=lambda x:history.append(objective(x)[0]),
                       options={"maxiter":max_iter, "ftol":min(tol*.01, 1e-12),
                                "gtol":min(tol, 1e-7), "maxls":40})
        grad = opt.jac.copy()
        for j, (lower, upper) in enumerate(bounds):
            if lower is not None and lower == upper:
                grad[j] = 0.
            if lower is not None and opt.x[j] <= lower+1e-9 and grad[j] > 0:
                grad[j] = 0.
        pg = float(np.max(np.abs(grad)))
        fits.append((opt, history, pg))
    accepted = [f for f in fits if f[0].success and f[2] < 1e-4]
    opt, history, pg = min(accepted or fits, key=lambda f:f[0].fun)
    mu, L = unpack(opt.x)
    S = L@L.T
    T = S[:-1, :-1]
    b = np.linalg.solve(T, S[:-1, -1])
    a = mu[-1]-b@mu[:-1]
    tau = float(L[-1, -1]**2)
    gain = np.linalg.solve(S+C, np.broadcast_to(S, C.shape)).transpose(0, 2, 1)
    latent = mu+np.einsum("nij,nj->ni", gain, z-mu)
    posterior = S-gain@S
    posterior = .5*(posterior+posterior.transpose(0, 2, 1))
    latent_boundary = bool(np.linalg.eigvalsh(T).min() < 1e-7)
    converged = bool(opt.success and pg < 1e-4 and not latent_boundary and np.isfinite(opt.fun))
    diagnostics = dict(converged=converged, objective="gaussian_structural_joint_density",
                       projected_gradient=pg, starts_succeeded=len(accepted),
                       latent_covariance_boundary=latent_boundary,
                       intrinsic_variance_boundary=bool(L[-1, -1] < 1e-6),
                       likelihood_kind="joint observed Gaussian density",
                       posterior_uncertainty="conditional on fitted parameters; parameter uncertainty excluded")
    return b, a, latent, posterior, tau, S, mu, opt, history, diagnostics
