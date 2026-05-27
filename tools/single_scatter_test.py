#!/usr/bin/env python3
"""
Single-scattering kernel comparison for igrmonty and CoportS.

The test fixes a local fluid-frame photon energy and electron temperature,
then samples one Compton scattering many times with:

  1. igrmonty-style direct Lorentz boosts;
  2. CoportS-style electron-frame tetrads.

Output values are dimensionless energy gains A = E_out / E_in in the local
fluid frame.
"""

import argparse
import csv
import math
import random
from statistics import mean, median


HPLANCK = 6.62607015e-27
ME = 9.1093837e-28
CL = 2.99792458e10
MEC2 = ME * CL * CL


def mj_fdist(gamma, thetae):
    if gamma <= 1.0:
        return 0.0
    return gamma * gamma * math.sqrt(gamma * gamma - 1.0) * math.exp(-gamma / thetae)


def find_mj_fmax(thetae):
    lo = math.log(max(1.000001, 0.01 * thetae))
    hi = math.log(max(100.0, 1000.0 * thetae))
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - gr * (hi - lo)
    d = lo + gr * (hi - lo)
    for _ in range(120):
        if mj_fdist(math.exp(c), thetae) < mj_fdist(math.exp(d), thetae):
            lo = c
            c = d
            d = lo + gr * (hi - lo)
        else:
            hi = d
            d = c
            c = hi - gr * (hi - lo)
    return mj_fdist(math.exp(0.5 * (lo + hi)), thetae)


def sample_gamma_igrmonty(thetae, fmax):
    """igrmonty src/compton.c sample_beta_distr_num for thermal EDF."""
    if thetae < 0.01:
        return 1.000001
    lge_min = math.log(max(1.0, 0.01 * thetae))
    lge_max = math.log(max(100.0, 1000.0 * thetae))
    while True:
        gamma = math.exp(lge_min + (lge_max - lge_min) * random.random())
        if mj_fdist(gamma, thetae) / fmax >= random.random():
            return gamma


def chi_square(n):
    return sum(random.gauss(0.0, 1.0) ** 2 for _ in range(n))


def sample_gamma_coports(thetae):
    """CoportS sample.cpp SampleMJD."""
    sqrt_pi = math.sqrt(math.pi)
    x1 = sqrt_pi / 4.0
    x2 = math.sqrt(thetae) / (2.0 * math.sqrt(2.0))
    x3 = 3.0 * sqrt_pi * thetae / 8.0
    x4 = thetae * math.sqrt(thetae / 2.0)
    norm = x1 + x2 + x3 + x4
    c1 = x1 / norm
    c2 = c1 + x2 / norm
    c3 = c2 + x3 / norm
    while True:
        u = random.random()
        if u < c1:
            n = 3
        elif u < c2:
            n = 4
        elif u < c3:
            n = 5
        else:
            n = 6
        y = math.sqrt(chi_square(n) / 2.0)
        accept = math.sqrt(1.0 + 0.5 * thetae * y * y) / (
            1.0 + y * math.sqrt(0.5 * thetae)
        )
        if accept > random.random():
            return 1.0 + thetae * y * y


def sample_mu(beta):
    if beta < 1.0e-14:
        return 2.0 * random.random() - 1.0
    x = random.random()
    det = 1.0 + 2.0 * beta + beta * beta - 4.0 * beta * x
    return max(-1.0, min(1.0, (1.0 - math.sqrt(max(0.0, det))) / beta))


def total_kn(a):
    if a < 1.0e-3:
        return 1.0 - 2.0 * a
    a2 = a * a
    return (3.0 / (4.0 * a2)) * (
        2.0
        + a2 * (1.0 + a) / ((1.0 + 2.0 * a) * (1.0 + 2.0 * a))
        + (a2 - 2.0 * a - 2.0) / (2.0 * a) * math.log(1.0 + 2.0 * a)
    )


def sample_electron(gamma_sampler, thetae, eps, fmax=None):
    while True:
        gamma = gamma_sampler(thetae, fmax) if fmax is not None else gamma_sampler(thetae)
        beta = math.sqrt(max(0.0, 1.0 - 1.0 / (gamma * gamma)))
        mu = sample_mu(beta)
        e_electron_frame = eps * gamma * (1.0 - beta * mu)
        if total_kn(e_electron_frame) > random.random():
            phi = 2.0 * math.pi * random.random()
            st = math.sqrt(max(0.0, 1.0 - mu * mu))
            gb = gamma * beta
            return [
                gamma,
                gb * st * math.cos(phi),
                gb * st * math.sin(phi),
                gb * mu,
            ]


def boost(v, u):
    """igrmonty src/compton.c boost in orthonormal coordinates."""
    gamma = u[0]
    vel = math.sqrt(abs(1.0 - 1.0 / (gamma * gamma)))
    denom = gamma * vel + 1.0e-40
    n = [u[1] / denom, u[2] / denom, u[3] / denom]
    gm1 = gamma - 1.0
    out = [0.0, 0.0, 0.0, 0.0]
    out[0] = u[0] * v[0] - u[1] * v[1] - u[2] * v[2] - u[3] * v[3]
    for i in range(3):
        out[i + 1] = -u[i + 1] * v[0]
        for j in range(3):
            out[i + 1] += ((1.0 if i == j else 0.0) + n[i] * n[j] * gm1) * v[j + 1]
    return out


def sample_thomson_mu():
    while True:
        x1 = 2.0 * random.random() - 1.0
        x2 = 0.75 * random.random()
        if x2 < 0.375 * (1.0 + x1 * x1):
            return x1


def klein_nishina_diff(a, ap):
    ch = 1.0 + 1.0 / a - 1.0 / ap
    return (a / ap + ap / a - 1.0 + ch * ch) / (a * a)


def sample_kn_energy_igrmonty(e):
    emin = e / (1.0 + 2.0 * e)
    emax = e
    ymax = 2.0 * (1.0 + 2.0 * e + 2.0 * e * e) / (e * e * (1.0 + 2.0 * e))
    while True:
        ep = emin + (emax - emin) * random.random()
        if ymax * random.random() < klein_nishina_diff(e, ep):
            return ep


def sample_kn_energy_coports(e):
    if e < 1.0:
        fmax = (2.0 + 4.0 * e + 4.0 * e * e) / (1.0 + 2.0 * e)
        lo = e / (1.0 + 2.0 * e)
        width = e - lo
        while True:
            ep = lo + width * random.random()
            f = e / ep + ep / e - 1.0 + (1.0 + 1.0 / e - 1.0 / ep) ** 2
            if fmax * random.random() < f:
                return ep
    lo = e / (1.0 + 2.0 * e)
    log_ratio = math.log(1.0 + 2.0 * e)
    while True:
        ep = lo * math.exp(random.random() * log_ratio)
        f = e / ep + ep / e - 1.0 + (1.0 + 1.0 / e - 1.0 / ep) ** 2
        if (2.0 * e / ep) * random.random() < f:
            return ep


def orthonormal_basis_from_z(v0):
    trial = [random.uniform(-1.0, 1.0) for _ in range(3)]
    dot = sum(trial[i] * v0[i] for i in range(3))
    v1 = [trial[i] - dot * v0[i] for i in range(3)]
    norm = math.sqrt(sum(x * x for x in v1))
    if norm < 1.0e-14:
        return orthonormal_basis_from_z(v0)
    v1 = [x / norm for x in v1]
    v2 = [
        v0[1] * v1[2] - v0[2] * v1[1],
        v0[2] * v1[0] - v0[0] * v1[2],
        v0[0] * v1[1] - v0[1] * v1[0],
    ]
    return v1, v2


def scatter_igrmonty(thetae, eps, fmax):
    k = [eps, 0.0, 0.0, eps]
    p = sample_electron(sample_gamma_igrmonty, thetae, eps, fmax)
    ke = boost(k, p)
    if ke[0] > 1.0e-4:
        ep = sample_kn_energy_igrmonty(ke[0])
        cth = 1.0 - 1.0 / ep + 1.0 / ke[0]
    else:
        ep = ke[0]
        cth = sample_thomson_mu()
    cth = max(-1.0, min(1.0, cth))
    sth = math.sqrt(max(0.0, 1.0 - cth * cth))
    kemag = math.sqrt(ke[1] * ke[1] + ke[2] * ke[2] + ke[3] * ke[3])
    v0 = [ke[1] / kemag, ke[2] / kemag, ke[3] / kemag]
    v1, v2 = orthonormal_basis_from_z(v0)
    phi = 2.0 * math.pi * random.random()
    direction = [
        cth * v0[i] + sth * (math.cos(phi) * v1[i] + math.sin(phi) * v2[i])
        for i in range(3)
    ]
    kpe = [ep, ep * direction[0], ep * direction[1], ep * direction[2]]
    p_inv = [p[0], -p[1], -p[2], -p[3]]
    kp = boost(kpe, p_inv)
    return kp[0] / eps


def mdot(a, b):
    return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]


def minkowski_orthogonal_tetrads(u, k, d):
    beta = mdot(d, u)
    omega = -mdot(k, u)
    cc = mdot(k, d) / omega - beta
    dnorm = mdot(d, d)
    nn2 = dnorm + beta * beta - cc * cc
    nn = math.sqrt(abs(nn2))
    e0 = u[:]
    e3 = [k[i] / omega - u[i] for i in range(4)]
    e2 = [(d[i] + beta * u[i] - cc * e3[i]) / nn for i in range(4)]
    v1 = d[3] * k[2] * u[1] - d[2] * k[3] * u[1] - d[3] * k[1] * u[2] + d[1] * k[3] * u[2] + d[2] * k[1] * u[3] - d[1] * k[2] * u[3]
    v2 = -d[3] * k[2] * u[0] + d[2] * k[3] * u[0] + d[3] * k[0] * u[2] - d[0] * k[3] * u[2] - d[2] * k[0] * u[3] + d[0] * k[2] * u[3]
    v3 = d[3] * k[1] * u[0] - d[1] * k[3] * u[0] - d[3] * k[0] * u[1] + d[0] * k[3] * u[1] + d[1] * k[0] * u[3] - d[0] * k[1] * u[3]
    v4 = -d[2] * k[1] * u[0] + d[1] * k[2] * u[0] + d[2] * k[0] * u[1] - d[0] * k[2] * u[1] - d[1] * k[0] * u[2] + d[0] * k[1] * u[2]
    scalar = 1.0 / (omega * nn)
    e1 = [scalar * (-v1), scalar * v2, scalar * v3, scalar * v4]
    return [e0, e1, e2, e3]


def scatter_coports(thetae, eps):
    pe = sample_electron(sample_gamma_coports, thetae, eps)
    k_local = [1.0, 0.0, 0.0, 1.0]
    e_ele = -eps * mdot(k_local, pe) / k_local[0]
    if e_ele < 1.0e-4:
        ep = e_ele
        cth = sample_thomson_mu()
    else:
        ep = sample_kn_energy_coports(e_ele)
        cth = 1.0 + 1.0 / e_ele - 1.0 / ep
    cth = max(-1.0, min(1.0, cth))
    sth = math.sqrt(max(0.0, 1.0 - cth * cth))
    phi = 2.0 * math.pi * random.random()
    k_es = [1.0, sth * math.cos(phi), sth * math.sin(phi), cth]
    tetrad = minkowski_orthogonal_tetrads(pe, k_local, [0.0, 1.0, 0.3, 0.01])
    k_us = [sum(tetrad[j][i] * k_es[j] for j in range(4)) for i in range(4)]
    e_us = ep * k_us[0] / k_es[0]
    return e_us / eps


def percentile(values, q):
    if not values:
        return float("nan")
    xs = sorted(values)
    pos = (len(xs) - 1) * q / 100.0
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    return xs[lo] * (hi - pos) + xs[hi] * (pos - lo)


def print_stats(name, values):
    print(f"{name}:")
    print(f"  n       = {len(values)}")
    print(f"  mean    = {mean(values):.8g}")
    print(f"  median  = {median(values):.8g}")
    for q in (90, 99, 99.9, 99.99):
        print(f"  p{q:<5} = {percentile(values, q):.8g}")
    print(f"  max     = {max(values):.8g}")
    for cut in (10, 100, 1000):
        frac = sum(1 for x in values if x > cut) / len(values)
        print(f"  P(A>{cut:<4}) = {frac:.8g}")


def maybe_plot(path, igr, co):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"plot skipped: matplotlib unavailable ({exc})")
        return
    bins = 120
    log_igr = [math.log10(x) for x in igr if x > 0.0]
    log_co = [math.log10(x) for x in co if x > 0.0]
    plt.figure(figsize=(8, 5))
    plt.hist(log_igr, bins=bins, histtype="step", density=True, label="igrmonty")
    plt.hist(log_co, bins=bins, histtype="step", density=True, label="CoportS")
    plt.xlabel(r"$\log_{10}(E_{\rm out}/E_{\rm in})$")
    plt.ylabel("PDF")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--thetae", type=float, default=16.0, help="dimensionless electron temperature")
    parser.add_argument("--eps", type=float, default=None, help="dimensionless photon energy h nu / mec^2")
    parser.add_argument("--nu", type=float, default=1.0e12, help="photon frequency in Hz, used if --eps is absent")
    parser.add_argument("--n", type=int, default=10000, help="number of single-scattering samples")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--out", default="single_scatter_test.csv")
    parser.add_argument("--plot", default=None, help="optional output PNG path")
    args = parser.parse_args()

    eps = args.eps if args.eps is not None else HPLANCK * args.nu / MEC2
    random.seed(args.seed)
    fmax = find_mj_fmax(args.thetae)

    igr = []
    co = []
    with open(args.out, "w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["sample", "gain_igrmonty", "gain_coports", "log10_gain_igrmonty", "log10_gain_coports"])
        for i in range(args.n):
            ai = scatter_igrmonty(args.thetae, eps, fmax)
            ac = scatter_coports(args.thetae, eps)
            igr.append(ai)
            co.append(ac)
            writer.writerow([i, ai, ac, math.log10(ai) if ai > 0 else "", math.log10(ac) if ac > 0 else ""])

    print(f"thetae = {args.thetae:g}")
    print(f"eps    = {eps:.8e}  (nu = {eps * MEC2 / HPLANCK:.8e} Hz)")
    print(f"wrote  {args.out}")
    print_stats("igrmonty", igr)
    print_stats("CoportS", co)
    print("ratio of p99 gains CoportS/igrmonty = %.8g" % (percentile(co, 99) / percentile(igr, 99)))
    print("ratio of p99.9 gains CoportS/igrmonty = %.8g" % (percentile(co, 99.9) / percentile(igr, 99.9)))
    if args.plot:
        maybe_plot(args.plot, igr, co)


if __name__ == "__main__":
    main()
