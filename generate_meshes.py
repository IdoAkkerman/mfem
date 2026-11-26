#!/usr/bin/env python3
"""
Generate 5 single-patch 2D NURBS meshes of the unit square for MFEM.

  square-nurbs-deformed-L0.mesh   eps = 0.0   (affine reference, straight square)
  square-nurbs-deformed-L1.mesh   eps = 0.2
  square-nurbs-deformed-L2.mesh   eps = 0.4
  square-nurbs-deformed-L3.mesh   eps = 0.6
  square-nurbs-deformed-L4.mesh   eps = 0.8

All levels: degree p = 2, 4x4 elements, uniform knots, unit weights.
Only the control points change between levels.

Analytic distortion (A = 2*pi*xi, B = 2*pi*eta):
    x = xi  + eps/(2*pi) * sin(A) * sin(B)
    y = eta + eps/(2*pi) * sin(A) * sin(B)
The bump is zero on the whole boundary, so the domain is always the unit
square. Exact det J = 1 + eps*sin(A+B) >= 1 - eps > 0 for eps < 1.

Usage:  python3 generate_meshes.py
"""
import numpy as np
from scipy.interpolate import BSpline

P, NEL = 2, 4
EPS = [0.0, 0.2, 0.4, 0.6, 0.8]

# open, uniform knot vector; n = NEL + P control points per direction
knots = np.r_[[0.0] * P, np.linspace(0, 1, NEL + 1), [1.0] * P]
n = len(knots) - P - 1
greville = np.array([knots[i + 1:i + P + 1].mean() for i in range(n)])


def F(xi, eta, eps):
    u = eps / (2 * np.pi) * np.sin(2 * np.pi * xi) * np.sin(2 * np.pi * eta)
    return xi + u, eta + u


def control_net(eps):
    X, Y = np.meshgrid(greville, greville)  # axis0 = eta, axis1 = xi
    return F(X, Y, eps)


def min_detJ(eps, ns=200):
    X, Y = control_net(eps)
    B = BSpline(knots, np.eye(n), P)
    s = np.linspace(0, 1 - 1e-12, ns)
    N, dN = B(s), B.derivative()(s)
    ev = lambda A, Bm, C: A @ C @ Bm.T
    xxi, xeta = ev(N, dN, X), ev(dN, N, X)
    yxi, yeta = ev(N, dN, Y), ev(dN, N, Y)
    return (xxi * yeta - xeta * yxi).min()


def write_mesh(path, eps):
    X, Y = control_net(eps)
    kv = " ".join(f"{k:g}" for k in knots)
    L = [
        "MFEM NURBS mesh v1.0", "",
        f"# unit square, p={P}, {NEL}x{NEL} elements, eps={eps}", "",
        "dimension", "2", "",
        "elements", "1", "1 3 0 1 2 3", "",
        "boundary", "4",
        "1 1 0 1", "2 1 1 2", "3 1 2 3", "4 1 3 0", "",
        "edges", "4",
        "0 0 1", "1 1 2", "0 3 2", "1 0 3", "",
        "vertices", "4", "",
        "patches", "",
        "# Patch 1", "",
        "knotvectors", "2",
        f"{P} {n} {kv}",
        f"{P} {n} {kv}", "",
        "dimension", "2", "",
        "controlpoints",
    ]
    for b in range(n):          # xi index runs fastest
        for a in range(n):
            L.append(f"{X[b, a]:.16g} {Y[b, a]:.16g} 1")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    for lvl, eps in enumerate(EPS):
        name = f"square-nurbs-deformed-L{lvl}.mesh"
        write_mesh(name, eps)
        print(f"{name}: eps={eps}, min detJ = {min_detJ(eps):.3f}")
