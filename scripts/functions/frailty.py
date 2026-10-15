"""Conditional gamma-frailty likelihood with a piecewise constant baseline.

Independent implementation. The review does not specify its baseline/profile
method, so these estimates must not be labelled exact reproductions.
"""

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.stats import chi2


def exposure_matrix(time, edges):
    return np.maximum(
        0, np.minimum(np.asarray(time)[:, None], np.asarray(edges)[1:]) - np.asarray(edges)[:-1]
    )


def gamma_loglik_gradient(parameters, x, entry_exposure, exit_exposure, event, event_band, theta):
    p = x.shape[1]
    beta = parameters[:p]
    lam = np.exp(parameters[p:])
    eta = x @ beta
    rr = np.exp(eta)
    he = entry_exposure @ lam
    hx = exit_exposure @ lam
    se = theta * rr * he
    sx = theta * rr * hx
    if theta == 0:
        survival = -rr * (hx - he)
    else:
        survival = (np.log1p(se) - np.log1p(sx)) / theta
    likelihood = np.sum(survival + event * (parameters[p:][event_band] + eta - np.log1p(sx)))
    ge = event / (1 + sx) - rr * hx / (1 + sx) + rr * he / (1 + se)
    coef_x = (1 + event * theta) * rr / (1 + sx)
    coef_e = rr / (1 + se)
    gb = (
        np.bincount(event_band, weights=event, minlength=len(lam))
        - (exit_exposure.T @ coef_x - entry_exposure.T @ coef_e) * lam
    )
    grad = np.r_[x.T @ ge, gb]
    return -float(likelihood), -grad


def fit_gamma(x, entry, exit, event, theta, edges=(18, 35, 45, 55, 65, 75, 100, 130), start=None):
    ae = exposure_matrix(entry, edges)
    ax = exposure_matrix(exit, edges)
    bands = np.clip(np.searchsorted(edges, exit, side='left') - 1, 0, len(edges) - 2)
    p = x.shape[1]
    if start is None:
        rates = np.bincount(bands, weights=event, minlength=len(edges) - 1) / np.maximum(
            (ax - ae).sum(axis=0), 1
        )
        start = np.r_[np.zeros(p), np.log(np.maximum(rates, 1e-7))]
    result = minimize(
        gamma_loglik_gradient,
        start,
        args=(x, ae, ax, np.asarray(event, float), bands, theta),
        jac=True,
        method='L-BFGS-B',
        bounds=[(-5, 5)] * p + [(-20, 0)] * (len(edges) - 1),
        options={'maxiter': 1000, 'ftol': 1e-12, 'gtol': 1e-5, 'maxls': 40},
    )
    if not result.success:
        raise RuntimeError('Gamma-frailty likelihood optimisation failed: ' + result.message)
    return result


def profile_theta(x, entry, exit, event):
    """Profile theta; boundary-aware interval reported as a modelling sensitivity."""
    fits = {}

    def at(theta):
        key = round(float(theta), 8)
        if key not in fits:
            start = None if not fits else min(fits.items(), key=lambda kv: abs(kv[0] - key))[1].x
            fits[key] = fit_gamma(x, entry, exit, event, key, start=start)
        return fits[key]

    from scipy.optimize import minimize_scalar

    at(0.0)
    optimum = minimize_scalar(
        lambda t: at(t).fun, bounds=(0, 3), method='bounded', options={'xatol': 0.002}
    )
    theta = float(optimum.x)
    if at(0).fun < at(theta).fun:
        theta = 0.0
    best = at(theta)
    cut = chi2.ppf(0.95, 1)

    def root(t):
        return 2 * (at(t).fun - best.fun) - cut

    lower = 0.0 if theta == 0 or root(0) < 0 else brentq(root, 0, theta, xtol=0.002)
    upper = None if root(3) < 0 else brentq(root, theta, 3, xtol=0.002)
    return theta, lower, upper, best, sorted((t, float(f.fun)) for t, f in fits.items())
