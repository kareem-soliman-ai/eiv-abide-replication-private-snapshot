"""Classical errors-in-variables comparison methods.

Implements OLS, WLS, Deming (constant error ratio), York (heteroscedastic
per-point weights, Williamson 1968 / York et al. 2004), SIMEX (Cook and
Stefanski 1994) and regression calibration. Shared interface: each `fit_*`
function returns (slope, intercept) or a small dataclass with convergence
diagnostics where relevant.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import linregress


@dataclass
class YorkFitResult:
    slope: float
    intercept: float
    iterations: int
    converged: bool


@dataclass
class SimexFitResult:
    slope: float
    intercept: float
    lambda_slopes: dict[float, float]
    quad_coeffs_slope: list[float]


def fit_ols(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    result = linregress(x, y)
    return float(result.slope), float(result.intercept)


def fit_wls(
    x: np.ndarray,
    y: np.ndarray,
    sigma_y: np.ndarray,
) -> tuple[float, float]:
    """Weighted least squares with w_i = 1 / sigma_y_i^2."""

    weights = 1.0 / (sigma_y**2)
    weight_sum = float(np.sum(weights))
    x_bar = float(np.sum(weights * x) / weight_sum)
    y_bar = float(np.sum(weights * y) / weight_sum)
    slope = float(
        np.sum(weights * (x - x_bar) * (y - y_bar))
        / np.sum(weights * (x - x_bar) ** 2)
    )
    intercept = float(y_bar - slope * x_bar)
    return slope, intercept


def fit_deming(
    x: np.ndarray,
    y: np.ndarray,
    delta: float,
) -> tuple[float, float]:
    """Deming regression with constant error variance ratio delta = var(eps_y)/var(eps_x).

    Closed-form solution (Deming 1943; see Linnet 1990 for the modern
    biostatistical treatment). Assumes homoscedastic errors.
    """

    x_mean = float(np.mean(x))
    y_mean = float(np.mean(y))
    x_centered = x - x_mean
    y_centered = y - y_mean

    s_xx = float(np.mean(x_centered**2))
    s_yy = float(np.mean(y_centered**2))
    s_xy = float(np.mean(x_centered * y_centered))

    if np.isclose(s_xy, 0.0):
        return 0.0, y_mean

    diff = s_yy - delta * s_xx
    slope = (diff + np.sqrt(diff**2 + 4.0 * delta * s_xy**2)) / (2.0 * s_xy)
    intercept = y_mean - slope * x_mean
    return float(slope), float(intercept)


def fit_york(
    x: np.ndarray,
    y: np.ndarray,
    sigma_x: np.ndarray,
    sigma_y: np.ndarray,
    tolerance: float = 1e-10,
    max_iterations: int = 500,
) -> YorkFitResult:
    """York et al. (2004) iterative heteroscedastic EIV regression.

    Implements the case of uncorrelated errors in X and Y with per-point
    sigma_x and sigma_y known. This is the unweighted-correlation form
    (r_i = 0) of the general York estimator.
    """

    w_x = 1.0 / (sigma_x**2)
    w_y = 1.0 / (sigma_y**2)

    slope, intercept = fit_ols(x, y)
    converged = False

    for iteration in range(1, max_iterations + 1):
        # York (2004), r_i=0: W_i = 1/(sy_i^2 + b^2 sx_i^2).
        weights = 1.0 / (sigma_y**2 + slope**2 * sigma_x**2)
        weight_sum = float(np.sum(weights))
        x_bar = float(np.sum(weights * x) / weight_sum)
        y_bar = float(np.sum(weights * y) / weight_sum)

        u = x - x_bar
        v = y - y_bar
        beta_i = weights * (u / w_y + slope * v / w_x)

        numerator = float(np.sum(weights * beta_i * v))
        denominator = float(np.sum(weights * beta_i * u))
        if np.isclose(denominator, 0.0):
            break
        new_slope = numerator / denominator

        if abs(new_slope - slope) < tolerance:
            slope = new_slope
            final_weights = 1.0 / (sigma_y**2 + slope**2 * sigma_x**2)
            intercept = float(np.sum(final_weights * (y - slope*x)) / np.sum(final_weights))
            converged = True
            return YorkFitResult(slope, intercept, iteration, True)

        slope = new_slope
        intercept = float(y_bar - slope * x_bar)

    return YorkFitResult(slope, intercept, iteration, converged)


def fit_simex(
    x: np.ndarray,
    y: np.ndarray,
    sigma_x: np.ndarray,
    n_sim: int = 200,
    lambdas: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0),
    seed: int = 42,
) -> SimexFitResult:
    """SIMEX (Cook and Stefanski 1994) with quadratic extrapolation to lambda=-1."""

    rng = np.random.default_rng(seed)
    n = len(x)

    ols_slope, ols_intercept = fit_ols(x, y)
    lambda_grid = (0.0, *lambdas)
    mean_slopes: dict[float, float] = {0.0: ols_slope}
    mean_intercepts: dict[float, float] = {0.0: ols_intercept}

    for lam in lambdas:
        slopes = np.empty(n_sim, dtype=float)
        intercepts = np.empty(n_sim, dtype=float)
        for i in range(n_sim):
            noise = rng.standard_normal(n) * sigma_x * np.sqrt(lam)
            x_star = x + noise
            s, c = fit_ols(x_star, y)
            slopes[i] = s
            intercepts[i] = c
        mean_slopes[float(lam)] = float(slopes.mean())
        mean_intercepts[float(lam)] = float(intercepts.mean())

    lam_arr = np.array(lambda_grid, dtype=float)
    slope_arr = np.array([mean_slopes[lam] for lam in lambda_grid], dtype=float)
    intercept_arr = np.array([mean_intercepts[lam] for lam in lambda_grid], dtype=float)

    coeffs_slope = np.polyfit(lam_arr, slope_arr, 2)
    coeffs_intercept = np.polyfit(lam_arr, intercept_arr, 2)
    corrected_slope = float(np.polyval(coeffs_slope, -1.0))
    corrected_intercept = float(np.polyval(coeffs_intercept, -1.0))

    return SimexFitResult(
        slope=corrected_slope,
        intercept=corrected_intercept,
        lambda_slopes=mean_slopes,
        quad_coeffs_slope=coeffs_slope.tolist(),
    )


def fit_regression_calibration(
    x: np.ndarray,
    y: np.ndarray,
    sigma_x: np.ndarray,
) -> tuple[float, float]:
    """Classical regression calibration using the reliability ratio.

    beta_RC = beta_OLS / kappa where kappa = var(X_true) / var(X_obs)
    and var(X_true) is estimated by var(X_obs) - mean(sigma_x^2).
    Assumes classical additive measurement error in X only.
    """

    ols_slope, _ = fit_ols(x, y)
    var_x_obs = float(np.var(x, ddof=1))
    mean_sigma_x_sq = float(np.mean(sigma_x**2))
    var_x_true = max(var_x_obs - mean_sigma_x_sq, 1e-12)
    kappa = var_x_true / var_x_obs
    slope_rc = ols_slope / kappa
    intercept_rc = float(np.mean(y) - slope_rc * np.mean(x))
    return float(slope_rc), intercept_rc
