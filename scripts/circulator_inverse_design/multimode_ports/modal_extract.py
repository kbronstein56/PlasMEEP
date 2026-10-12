"""Power-normalized modal a,b extraction on a port line (matches modal_s_core signs)."""
from __future__ import annotations

import numpy as np

from modes import normalize_on_line, power_scale, transverse_mode


def sample_hz_dn(sampler, hz, tris, xy, n_hat):
    n_hat = np.asarray(n_hat, float)
    n_hat = n_hat / np.linalg.norm(n_hat)
    xy = np.atleast_2d(np.asarray(xy, float))
    tid = sampler.locate(xy)
    Hz = np.zeros(len(xy), dtype=np.complex128)
    dn = np.zeros(len(xy), dtype=np.complex128)
    for i in range(len(xy)):
        t = int(tid[i])
        if t < 0:
            j = int(np.argmin(np.sum((sampler.points - xy[i]) ** 2, axis=1)))
            Hz[i] = hz[j]
            continue
        nodes = tris[t]
        ell = sampler._bary(t, xy[i])
        un = hz[nodes]
        Hz[i] = np.dot(ell, un)
        gx = np.dot(sampler.bx[t], un)
        gy = np.dot(sampler.by[t], un)
        dn[i] = n_hat[0] * gx + n_hat[1] * gy
    return Hz, dn


def extract_mode_m(
    sampler,
    hz,
    tris,
    line: dict,
    k0: float,
    m: int,
    W: float,
) -> tuple[complex, complex]:
    """Return power-normalized (a, b) for transverse mode m."""
    weights = np.asarray(line["weights"], float)
    if "offsets" in line:
        off = np.asarray(line["offsets"], float)
    else:
        span = float(np.sum(weights))
        off = np.linspace(-span / 2.0, span / 2.0, len(weights))
    tang = np.asarray(line["tangent"], float)
    tang = tang / np.linalg.norm(tang)
    mon = np.asarray(line["monitor_xy"], float)
    # y along aperture = tang · (xy - center)
    center = mon[len(mon) // 2]
    y = np.array([np.dot(p - center, tang) for p in mon])
    phi, _ = normalize_on_line(y, weights, m, W)
    n_hat = np.asarray(line["outward"], float)
    Hz, dn = sample_hz_dn(sampler, hz, tris, mon, n_hat)
    disc = float(k0**2 - (m * np.pi / W) ** 2)
    if disc > 0:
        beta = complex(np.sqrt(disc))
    else:
        alpha = float(np.sqrt(max(-disc, 0.0)))
        beta = 1j * max(alpha, 1e-30)
    c_H = complex(np.sum(weights * phi * Hz))
    c_d = complex(np.sum(weights * phi * dn))
    denom = 1j * beta
    a_f = 0.5 * (c_H - c_d / denom)
    b_f = 0.5 * (c_H + c_d / denom)
    sc = power_scale(beta, k0)
    return a_f * sc, b_f * sc
