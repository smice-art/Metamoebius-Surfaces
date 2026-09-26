bl_info = {
    "name": "Metamobius Surfaces (Checkerboard Surface Generator)",
    "author": "Claudio Wabner & Claude",
    "version": (2, 3, 0),
    "blender": (4, 5, 0),
    "location": "View3D > N-Panel > Metamobius",
    "description": "One-sided (non-orientable) spanning surfaces built from symmetric "
                   "knot diagrams: shaded regions become membranes, every crossing a "
                   "half-twisted band (checkerboard / Tait surface construction)",
    "category": "Add Mesh",
}

# =============================================================================
#  HOW IT WORKS
# -----------------------------------------------------------------------------
#  1. A closed planar curve (or several) is drawn -> the knot DIAGRAM.
#  2. The regions of the diagram are 2-coloured like a chessboard. The colour of
#     a point is the parity of the winding number (even-odd rule).
#  3. The regions of one colour are triangulated -> the MEMBRANES.
#  4. A small disk is cut out around every crossing; inside it a HALF-TWISTED
#     BAND connects the two shaded corners (over strand at +h, under at -h).
#  5. Membrane heights are solved as a harmonic / Poisson field so everything
#     is smooth; an optional 3D relax gives a soap-film look.
#  The result is a spanning surface of the knot. For alternating diagrams one of
#  the two colourings is almost always non-orientable (one-sided).
# =============================================================================

import math
import numpy as np
import bpy
from mathutils import Vector
from mathutils.geometry import delaunay_2d_cdt
from bpy.props import (EnumProperty, IntProperty, FloatProperty, BoolProperty,
                       PointerProperty, StringProperty)
from bpy.types import PropertyGroup, Operator, Panel

OBJ_NAME = "Metamobius_Surface"
EDGE_NAME = "Metamobius_Boundary"
MAT_NAME = "Metamobius_Material"
EDGE_MAT_NAME = "Metamobius_Edge_Material"


# =============================================================================
#  1. DIAGRAM CURVES
# =============================================================================

def _gcd(*v):
    g = 0
    for x in v:
        g = math.gcd(g, abs(int(x)))
    return max(g, 1)


def diagram_components(P):
    """Returns a list of closed polylines (N,2) - dense, not yet resampled."""
    fam = P["family"]
    n = 6000
    comps = []
    if fam == 'ROSETTE':
        W, L = P["windings"], P["lobes"]
        g = _gcd(W, L)
        W, L = W // g, L // g
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        r = 1.0 + P["amp"] * np.cos(L * t)
        comps.append(np.c_[r * np.cos(W * t), r * np.sin(W * t)])
    elif fam == 'EPI':
        f = [P["f1"], P["f2"], P["f3"]]
        a = [1.0, P["a2"], P["a3"]]
        ph = [0.0, 0.0, P["ph3"]]      # with 3 terms only ONE phase changes the shape
        g = _gcd(*[fi for fi, ai in zip(f, a) if abs(ai) > 1e-9])
        t = np.linspace(0, 2 * np.pi / g, n, endpoint=False)
        z = sum(ai * np.exp(1j * (fi * t + pi)) for fi, ai, pi in zip(f, a, ph))
        comps.append(np.c_[z.real, z.imag])
    elif fam == 'LISSAJOUS':
        a, b = P["lx"], P["ly"]
        g = _gcd(a, b)
        a, b = a // g, b // g
        t = np.linspace(0, 2 * np.pi, n, endpoint=False)
        comps.append(np.c_[np.sin(a * t + P["lphase"]), P["laspect"] * np.sin(b * t)])
    elif fam == 'RINGS':
        m = P["rings"]
        t = np.linspace(0, 2 * np.pi, max(800, n // m), endpoint=False)
        for k in range(m):
            ang = 2 * np.pi * k / m
            cx, cy = P["ring_dist"] * np.cos(ang), P["ring_dist"] * np.sin(ang)
            rx = P["ring_radius"]
            ry = P["ring_radius"] * P["ring_squash"]
            x, y = rx * np.cos(t), ry * np.sin(t)
            ca, sa = np.cos(ang), np.sin(ang)
            comps.append(np.c_[cx + ca * x - sa * y, cy + sa * x + ca * y])
    elif fam == 'CHAIN':
        m = P["rings"]
        t = np.linspace(0, 2 * np.pi, max(800, n // m), endpoint=False)
        s = P["chain_spacing"]
        for k in range(m):
            cx = (k - (m - 1) / 2.0) * s
            comps.append(np.c_[cx + np.cos(t), P["ring_squash"] * np.sin(t)])
    return comps


def normalize_components(comps, size=1.0):
    allp = np.vstack(comps)
    c = 0.5 * (allp.min(0) + allp.max(0))
    rad = np.max(np.linalg.norm(allp - c, axis=1))
    return [(p - c) / rad * size for p in comps]


def resample(poly, ds):
    seg = np.roll(poly, -1, axis=0) - poly
    ln = np.linalg.norm(seg, axis=1)
    cum = np.r_[0, np.cumsum(ln)]
    total = cum[-1]
    m = max(24, int(total / ds))
    s = np.linspace(0, total, m, endpoint=False)
    idx = np.searchsorted(cum, s, side='right') - 1
    idx = np.clip(idx, 0, len(poly) - 1)
    u = (s - cum[idx]) / np.maximum(ln[idx], 1e-15)
    return poly[idx] + seg[idx] * u[:, None]


def invert_components(comps, p, r0=1.0):
    out = []
    for c in comps:
        d = c - p
        d2 = np.sum(d * d, axis=1, keepdims=True)
        out.append(p + r0 * r0 * d / np.maximum(d2, 1e-12))
    return out


# =============================================================================
#  2. GEOMETRY HELPERS (vectorised, chunked)
# =============================================================================

def segments_of(comps):
    A, B = [], []
    for c in comps:
        A.append(c)
        B.append(np.roll(c, -1, axis=0))
    return np.vstack(A), np.vstack(B)


def parity(points, SA, SB, chunk=4096):
    """Even-odd rule: 1 = odd winding (shaded), 0 = even."""
    points = np.atleast_2d(points)
    res = np.zeros(len(points), dtype=np.int64)
    y1, y2 = SA[:, 1], SB[:, 1]
    x1, x2 = SA[:, 0], SB[:, 0]
    dy = np.where(np.abs(y2 - y1) < 1e-300, 1e-300, y2 - y1)
    for i in range(0, len(points), chunk):
        P = points[i:i + chunk]
        px, py = P[:, 0:1], P[:, 1:2]
        cond = (y1 > py) != (y2 > py)
        xi = x1 + (py - y1) * (x2 - x1) / dy
        res[i:i + chunk] = np.sum(cond & (px < xi), axis=1) & 1
    return res


def min_dist_to_points(points, cloud, chunk=2048):
    out = np.empty(len(points))
    for i in range(0, len(points), chunk):
        d = points[i:i + chunk, None, :] - cloud[None, :, :]
        out[i:i + chunk] = np.sqrt(np.min(np.sum(d * d, axis=2), axis=1))
    return out


def find_crossings(comps):
    """Returns list of dicts: pos, strands [(comp, seg, u), (comp, seg, u)]."""
    segs = []
    for ci, c in enumerate(comps):
        n = len(c)
        for si in range(n):
            segs.append((ci, si, n))
    segs = np.array(segs)
    SA, SB = segments_of(comps)
    R = SB - SA
    ns = len(SA)
    out = []
    chunk = 256
    for i0 in range(0, ns, chunk):
        i1 = min(ns, i0 + chunk)
        p, r = SA[i0:i1, None, :], R[i0:i1, None, :]
        q, s = SA[None, :, :], R[None, :, :]
        rxs = r[..., 0] * s[..., 1] - r[..., 1] * s[..., 0]
        qp = q - p
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (qp[..., 0] * s[..., 1] - qp[..., 1] * s[..., 0]) / rxs
            u = (qp[..., 0] * r[..., 1] - qp[..., 1] * r[..., 0]) / rxs
        hit = (np.abs(rxs) > 1e-14) & (t >= 0) & (t < 1) & (u >= 0) & (u < 1)
        ii, jj = np.nonzero(hit)
        for a, b in zip(ii, jj):
            gi, gj = i0 + a, b
            if gj <= gi:
                continue
            ci, si, n = segs[gi]
            cj, sj, _ = segs[gj]
            if ci == cj:
                d = abs(si - sj)
                if d <= 1 or d >= n - 1:
                    continue
            tt = t[a, b]
            out.append({"pos": SA[gi] + R[gi] * tt,
                        "strands": [(int(ci), int(si), float(tt)),
                                    (int(cj), int(sj), float(u[a, b]))]})
    return out


def exit_point(poly, seg, u, c, rho, direction):
    """Walk along the polyline from (seg,u) until distance to c == rho.
    Returns (point, arc-index-float, list of sample indices inside the disk)."""
    n = len(poly)
    inside = []
    if direction > 0:
        prev = poly[seg] + (poly[(seg + 1) % n] - poly[seg]) * u
        j = (seg + 1) % n
        for _ in range(n):
            if np.linalg.norm(poly[j] - c) >= rho:
                break
            inside.append(j)
            prev = poly[j]
            j = (j + 1) % n
        a, b = prev, poly[j]
        key = j - 0.5
    else:
        prev = poly[seg] + (poly[(seg + 1) % n] - poly[seg]) * u
        j = seg
        for _ in range(n):
            if np.linalg.norm(poly[j] - c) >= rho:
                break
            inside.append(j)
            prev = poly[j]
            j = (j - 1) % n
        a, b = prev, poly[j]
        key = j + 0.5
    # solve |a + s(b-a) - c| = rho
    d = b - a
    f = a - c
    A = d @ d
    B = 2 * f @ d
    C = f @ f - rho * rho
    disc = max(B * B - 4 * A * C, 0.0)
    s = (-B + math.sqrt(disc)) / (2 * A) if A > 1e-18 else 0.0
    s = min(max(s, 0.0), 1.0)
    return a + d * s, key % n, inside


def ang_norm(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


# =============================================================================
#  3. SURFACE BUILDER
# =============================================================================

def build_surface(P, report=print):
    size = 1.0
    g = 2.0 * size / max(8, P["density"])      # interior spacing (canvas units)
    ds = g * 0.45                              # boundary sample spacing
    sphere = P["canvas"] == 'SPHERE'

    comps = normalize_components(diagram_components(P), size)
    if abs(P["rotation"]) > 1e-9:
        c_, s_ = math.cos(P["rotation"]), math.sin(P["rotation"])
        Rm = np.array([[c_, s_], [-s_, c_]])
        comps = [c @ Rm for c in comps]
    comps = [resample(c, ds * 0.3) for c in comps]

    tp = 1                                      # shaded parity
    if P["shade"] == 'B':
        if not sphere:
            SA, SB = segments_of(comps)
            o = np.array([[0.0, 0.0]])
            if parity(o, SA, SB)[0] == 1 and min_dist_to_points(o, np.vstack(comps))[0] > 0.08:
                # invert about the centre: symmetric, outer region becomes bounded
                comps = invert_components(comps, o[0], 1.0)
                comps = normalize_components(comps, size)
                comps = [resample(c, ds * 0.3) for c in comps]
            else:
                sphere = True                   # outer region needs the sphere
                report("Colour B needs the sphere here - using Sphere canvas")
        if sphere:
            tp = 0                              # outer region = polar cap

    # ---- metric: plane length per canvas length -------------------------------
    s_st = 1.0 / math.tan(0.5 * math.pi * P["coverage"]) if sphere else 1.0

    def lam(p):
        p = np.atleast_2d(p)
        if not sphere:
            return np.ones(len(p))
        r2 = np.sum(p * p, axis=1) / (s_st * s_st)
        return s_st * (1.0 + r2) * 0.5

    def resample_metric(poly, step):
        seg = np.roll(poly, -1, axis=0) - poly
        ln = np.linalg.norm(seg, axis=1) / lam(poly + 0.5 * seg)
        cum = np.r_[0, np.cumsum(ln)]
        m = max(24, int(cum[-1] / step))
        sgrid = np.linspace(0, cum[-1], m, endpoint=False)
        idx = np.clip(np.searchsorted(cum, sgrid, side='right') - 1, 0, len(poly) - 1)
        u = (sgrid - cum[idx]) / np.maximum(ln[idx], 1e-15)
        return poly[idx] + seg[idx] * u[:, None]

    comps = [resample_metric(c, ds) for c in comps]
    SA, SB = segments_of(comps)
    allpts = np.vstack(comps)
    allpts_lam = lam(allpts)

    # ---- crossings -------------------------------------------------------------
    X = find_crossings(comps)
    if not X:
        raise RuntimeError("This diagram has no crossings - change the parameters.")
    cpos = np.array([x["pos"] for x in X])
    rho0 = P["band_size"] * 0.12 * size
    cumlen = []
    for c in comps:
        ln = np.linalg.norm(np.roll(c, -1, 0) - c, axis=1)
        cumlen.append(np.r_[0, np.cumsum(ln)])

    for k, x in enumerate(X):
        c = x["pos"]
        lc = lam(c)[0]
        dc = np.linalg.norm(cpos - c, axis=1)
        dc[k] = np.inf
        rho = min(rho0 * lc, 0.42 * dc.min())
        mask = np.ones(len(allpts), bool)
        off = 0
        for ci, cc in enumerate(comps):
            n = len(cc)
            for (sc, si, u) in x["strands"]:
                if sc != ci:
                    continue
                s0 = cumlen[ci][si]
                tot = cumlen[ci][-1]
                dd = np.abs(cumlen[ci][:-1] - s0)
                dd = np.minimum(dd, tot - dd)
                mask[off:off + n] &= dd > 3.0 * rho0 * lc
            off += n
        if mask.any():
            rho = min(rho, 0.7 * np.min(np.linalg.norm(allpts[mask] - c, axis=1)))
        x["rho"] = max(rho, ds * lc * 1.2)

    # ---- exits & shaded sectors ---------------------------------------------
    removed = [set() for _ in comps]
    bad = 0
    for x in X:
        c, rho = x["pos"], x["rho"]
        ex = []
        for sidx, (ci, si, u) in enumerate(x["strands"]):
            for d in (+1, -1):
                p, key, inside = exit_point(comps[ci], si, u, c, rho, d)
                removed[ci].update(inside)
                ex.append({"strand": sidx, "dir": d, "p": p, "key": key, "comp": ci,
                           "ang": math.atan2(p[1] - c[1], p[0] - c[0])})
        ex.sort(key=lambda e: e["ang"])
        x["exits"] = ex
        if ex[0]["strand"] == ex[1]["strand"] or ex[1]["strand"] == ex[2]["strand"]:
            bad += 1
        mids = []
        for i in range(4):
            a0 = ex[i]["ang"]
            da = (ex[(i + 1) % 4]["ang"] - a0) % (2 * np.pi)
            am = a0 + 0.5 * da
            mids.append(c + 0.5 * rho * np.array([math.cos(am), math.sin(am)]))
        par = (parity(np.array(mids), SA, SB) == tp).astype(int)
        x["sectors"] = (0, 2) if par[0] + par[2] >= par[1] + par[3] else (1, 3)
    if bad:
        raise RuntimeError(f"{bad} crossing(s) too tangent - adjust parameters or Band Size.")

    # ---- over / under ------------------------------------------------------
    rule = P["crossing_rule"]
    rng = np.random.default_rng(P["seed"])
    for x in X:
        ex = x["exits"]
        s_start = ex[x["sectors"][0]]["strand"]
        over = s_start if not P["handedness"] else 1 - s_start
        if rule == 'TORUS' and P["family"] == 'ROSETTE':
            vals = []
            for (ci, si, u) in x["strands"]:
                cc = comps[ci]
                d = cc[(si + 1) % len(cc)] - cc[si - 1]
                p = cc[si]
                vals.append(-(p @ d) / max(np.linalg.norm(p), 1e-9))  # ~ sin(L t)
            over = 0 if vals[0] > vals[1] else 1
            if P["handedness"]:
                over = 1 - over
        elif rule == 'RANDOM':
            over = int(rng.integers(0, 2))
        x["over"] = over

    # ---- vertices / constraints ---------------------------------------------
    h = P["twist_height"] * size * 0.12
    verts2, vz, vkind, vdisk, cons = [], [], [], [], []

    def addv(p, z, kind, disk=-1):
        verts2.append((float(p[0]), float(p[1])))
        vz.append(z); vkind.append(kind); vdisk.append(disk)
        return len(verts2) - 1

    ev = [[] for _ in comps]
    for xi, x in enumerate(X):
        for e in x["exits"]:
            e["z"] = (1.0 if e["strand"] == x["over"] else -1.0) * h
            e["vid"] = addv(e["p"], e["z"], 0, xi)
            ev[e["comp"]].append((e["key"], e))

    for ci, cc in enumerate(comps):
        n = len(cc)
        items = [(float(i), "s", i) for i in range(n) if i not in removed[ci]]
        items += [(key, "e", e) for (key, e) in ev[ci]]
        items.sort(key=lambda it: it[0])
        if not ev[ci]:
            ids = [addv(cc[i], 0.0, 0) for (_, _, i) in items]
            for a in range(len(ids)):
                cons.append((ids[a], ids[(a + 1) % len(ids)]))
            continue
        start = next(k for k, it in enumerate(items) if it[1] == "e" and it[2]["dir"] > 0)
        items = items[start:] + items[:start]
        chains, cur = [], []
        for it in items:
            cur.append(it)
            if it[1] == "e" and it[2]["dir"] < 0:
                chains.append(cur)
                cur = []
        for ch in chains:
            z0, z1 = ch[0][2]["z"], ch[-1][2]["z"]
            pts = np.array([cc[it[2]] if it[1] == "s" else it[2]["p"] for it in ch])
            L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
            tot = max(L[-1], 1e-12)
            ids = []
            for k, it in enumerate(ch):
                if it[1] == "e":
                    ids.append(it[2]["vid"])
                else:
                    s = L[k] / tot
                    sm = s * s * (3 - 2 * s)
                    zz = z0 * (1 - sm) + z1 * sm
                    zz *= 1.0 - P["strand_flatten"] * math.sin(math.pi * s)
                    ids.append(addv(pts[k], zz, 0))
            for a in range(len(ids) - 1):
                cons.append((ids[a], ids[a + 1]))

    # ---- band arcs ------------------------------------------------------------
    K = max(4, int(P["band_res"]))
    bands = []
    for xi, x in enumerate(X):
        c, rho, ex = x["pos"], x["rho"], x["exits"]
        arcs = []
        for si in x["sectors"]:
            e_a, e_b = ex[si], ex[(si + 1) % 4]
            a0 = e_a["ang"]
            da = (e_b["ang"] - a0) % (2 * np.pi)
            if e_a["strand"] == x["over"]:
                ang = [a0 + da * k / K for k in range(K + 1)]
                e_first, e_last = e_a, e_b
            else:
                ang = [a0 + da * (1 - k / K) for k in range(K + 1)]
                e_first, e_last = e_b, e_a
            ids = [e_first["vid"]]
            for k in range(1, K):
                p = c + rho * np.array([math.cos(ang[k]), math.sin(ang[k])])
                ids.append(addv(p, h * math.cos(math.pi * k / K), 1, xi))
            ids.append(e_last["vid"])
            for k in range(K):
                cons.append((ids[k], ids[k + 1]))
            arcs.append(ids)
        bands.append(arcs)

    # ---- polar cap (sphere, outer region shaded) ------------------------------
    cap = sphere and tp == 0
    Rcap = None
    if cap:
        eps = 1.5 * g
        Rcap = s_st * math.tan(0.5 * (math.pi - eps))
        ncap = max(12, int(2 * math.pi * math.sin(eps) / (0.8 * g)))
        cid = [addv((Rcap * math.cos(2 * math.pi * k / ncap), Rcap * math.sin(2 * math.pi * k / ncap)),
                    np.nan, 2) for k in range(ncap)]
        for k in range(ncap):
            cons.append((cid[k], cid[(k + 1) % ncap]))

    # ---- interior points ----------------------------------------------------
    if sphere:
        npts = int(4 * math.pi / (0.866 * g * g))
        i = np.arange(npts) + 0.5
        th = np.arccos(1 - 2 * i / npts)            # from south pole (origin)
        ph = math.pi * (1 + 5 ** 0.5) * i
        keep = th < math.pi - 2.2 * g
        th, ph = th[keep], ph[keep]
        rr = s_st * np.tan(0.5 * th)
        grid = np.c_[rr * np.cos(ph), rr * np.sin(ph)]
    else:
        xs = np.arange(-size - g, size + g, g)
        GX, GY = np.meshgrid(xs, xs * 0.866)
        GX = GX + (np.arange(GX.shape[0])[:, None] % 2) * 0.5 * g
        grid = np.c_[GX.ravel(), GY.ravel()]
    grid = grid[parity(grid, SA, SB) == tp]
    if len(grid):
        lg = lam(grid)
        grid, lg = grid[min_dist_to_points(grid, allpts) / lg > 0.55 * g], lg[min_dist_to_points(grid, allpts) / lg > 0.55 * g]
    if len(grid):
        rhos = np.array([x["rho"] for x in X])
        dx = np.min(np.linalg.norm(grid[:, None, :] - cpos[None, :, :], axis=2) - rhos[None, :], axis=1)
        grid = grid[dx / lg > 0.5 * g]
    for p in grid:
        addv(p, np.nan, 2)

    # ---- triangulate -----------------------------------------------------------
    vco = [Vector(v) for v in verts2]
    out = delaunay_2d_cdt(vco, cons, [], 0, 1e-10, True)
    ov, of, orig_v = out[0], out[2], out[3]
    nin = len(verts2)
    in2out = np.full(nin, -1, np.int64)
    for oi, lst in enumerate(orig_v):
        for ii in lst:
            if ii < nin:
                in2out[ii] = oi
    if (in2out < 0).any():
        raise RuntimeError("Triangulation merged vertices - lower Mesh Density or change Band Size.")
    OV = np.array([(v.x, v.y) for v in ov])
    OZ = np.full(len(OV), np.nan)
    OK = np.full(len(OV), 2, np.int64)
    OD = np.full(len(OV), -1, np.int64)
    for ii in range(nin):
        o = in2out[ii]
        if not np.isnan(vz[ii]):
            OZ[o] = vz[ii]
        OK[o] = min(OK[o], vkind[ii])
        if vdisk[ii] >= 0:
            OD[o] = vdisk[ii]
    faces = np.array([f for f in of if len(f) == 3], dtype=np.int64)
    d3 = OD[faces]
    faces = faces[~((d3[:, 0] >= 0) & (d3[:, 0] == d3[:, 1]) & (d3[:, 1] == d3[:, 2]))]
    cen = OV[faces].mean(axis=1)
    keep = parity(cen, SA, SB) == tp
    if cap:
        keep &= np.linalg.norm(cen, axis=1) < Rcap
    faces = faces[keep]

    V2 = [p for p in OV]
    Z = list(OZ)
    KIND = list(OK)
    extra = []
    extra_lab = []
    if cap:
        pole = len(V2)
        V2.append(np.array([np.inf, np.inf])); Z.append(np.nan); KIND.append(2)
        cc = [in2out[i] for i in cid]
        for k in range(len(cc)):
            extra.append((cc[k], cc[(k + 1) % len(cc)], pole))
            extra_lab.append(-1)

    # ---- twisted band patches ------------------------------------------------
    M = max(4, int(P["band_res"]))
    for xi, arcs in enumerate(bands):
        A = [in2out[i] for i in arcs[0]]
        B = [in2out[i] for i in arcs[1]]
        rows = []
        for k in range(K + 1):
            pa, pb = OV[A[k]], OV[B[k]]
            zk = h * math.cos(math.pi * k / K)
            row = [A[k]]
            for m in range(1, M):
                V2.append(pa + (pb - pa) * (m / M)); Z.append(zk)
                KIND.append(0 if k in (0, K) else 3)
                row.append(len(V2) - 1)
            row.append(B[k])
            rows.append(row)
        for k in range(K):
            r0, r1 = rows[k], rows[k + 1]
            for m in range(M):
                extra.append((r0[m], r0[m + 1], r1[m + 1]))
                extra.append((r0[m], r1[m + 1], r1[m]))
                extra_lab += [xi, xi]

    V2 = np.array(V2)
    Z = np.array(Z, dtype=float)
    KIND = np.array(KIND)
    F = np.vstack([faces, np.array(extra, dtype=np.int64)]) if extra else faces
    FL = np.r_[np.full(len(faces), -1, np.int64), np.array(extra_lab, dtype=np.int64)]
    used = np.zeros(len(V2), bool)
    used[F.ravel()] = True
    remap = -np.ones(len(V2), np.int64)
    remap[used] = np.arange(used.sum())
    V2, Z, KIND, F = V2[used], Z[used], KIND[used], remap[F]

    # ---- canvas embedding --------------------------------------------------
    if sphere:
        fin = np.isfinite(V2[:, 0])
        q = np.where(fin[:, None], V2, 0.0) / s_st
        r2 = np.sum(q * q, axis=1)
        # rotated by pi about X so the diagram centre faces +Z (the viewer)
        S = np.c_[2 * q[:, 0] / (1 + r2), -2 * q[:, 1] / (1 + r2), (1 - r2) / (1 + r2)]
        S[~fin] = (0.0, 0.0, -1.0)
        BASE, NRM = S, S
    else:
        BASE = np.c_[V2, np.zeros(len(V2))]
        NRM = np.tile([0.0, 0.0, 1.0], (len(V2), 1))

    # ---- heights: Poisson solve --------------------------------------------
    E = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    E = np.unique(np.sort(E, axis=1), axis=0)
    nv = len(V2)
    free = np.isnan(Z)
    deg = np.bincount(E.ravel(), minlength=nv).astype(float)

    def lap(x):
        s = np.bincount(E[:, 0], weights=x[E[:, 1]], minlength=nv) + \
            np.bincount(E[:, 1], weights=x[E[:, 0]], minlength=nv)
        return deg * x - s

    zfix = np.where(free, 0.0, Z)
    b = -lap(zfix)
    b[free] += P["inflate"] * 0.8 * g * g
    b[~free] = 0.0

    def Aop(x):
        y = np.where(free, x, 0.0)
        return np.where(free, lap(y), 0.0)

    x = np.zeros(nv); r = b - Aop(x); p = r.copy(); rs = r @ r
    for _ in range(4000):
        Ap = Aop(p)
        den = p @ Ap
        if abs(den) < 1e-30:
            break
        al = rs / den
        x += al * p; r -= al * Ap
        rn = r @ r
        if rn < 1e-22:
            break
        p = r + (rn / rs) * p; rs = rn
    Z = np.where(free, x, zfix)
    V = BASE + NRM * Z[:, None]

    # ---- boundary (knot) smoothing ------------------------------------------
    loops0 = boundary_loops(F)
    onb = np.zeros(nv, bool)
    for lp in loops0:
        onb[lp] = True
    for _ in range(int(P["edge_smooth"])):
        for lp in loops0:
            lp = np.array(lp)
            Q = V[lp]
            avg = 0.5 * (np.roll(Q, 1, 0) + np.roll(Q, -1, 0))
            V[lp] = Q + 0.5 * (avg - Q)

    # ---- light seam relax (smooths band/membrane junctions) -------------------
    mov = ~onb
    for _ in range(15):
        S3 = np.zeros_like(V)
        np.add.at(S3, E[:, 0], V[E[:, 1]])
        np.add.at(S3, E[:, 1], V[E[:, 0]])
        avg = S3 / np.maximum(deg, 1)[:, None]
        V[mov] += 0.5 * (avg[mov] - V[mov])

    # ---- soap film: blend towards the discrete minimal surface ---------------
    sf = float(P["soap_film"])
    if sf > 1e-6:
        Vmin = minimal_surface(V, F, onb, outer=int(P["soap_iter"]))
        V = V + sf * (Vmin - V)

    if not sphere and abs(P["dome"]) > 1e-9:
        V[:, 2] += P["dome"] * 0.5 * (1.0 - V[:, 0] ** 2 - V[:, 1] ** 2)

    V *= P["size"]
    F, orientable = orient_structured(F, FL)
    chi = nv - len(E) + len(F)
    loops = boundary_loops(F)
    if len(loops) > len(comps):
        raise RuntimeError(f"Surface broke up ({len(loops)} boundary loops for {len(comps)} knot "
                           "component(s)) - lower Band Size or raise Mesh Density.")
    info = {"crossings": len(X), "chi": int(chi), "boundaries": len(loops),
            "orientable": orientable, "verts": nv, "faces": len(F)}
    return V, F, loops, info




def cotan_edges(V, F):
    """Unique edges (E,2) and symmetric cotangent weights (clamped positive)."""
    i, j, k = F[:, 0], F[:, 1], F[:, 2]

    def cot(a, b, c):                     # cot of the angle at a
        u, v = V[b] - V[a], V[c] - V[a]
        cr = np.linalg.norm(np.cross(u, v), axis=1)
        return np.einsum('ij,ij->i', u, v) / np.maximum(cr, 1e-12)

    ca, cb, cc = cot(i, j, k), cot(j, k, i), cot(k, i, j)
    Eall = np.vstack([np.c_[j, k], np.c_[k, i], np.c_[i, j]])
    Wall = 0.5 * np.r_[ca, cb, cc]
    Eall = np.sort(Eall, axis=1)
    E, inv = np.unique(Eall, axis=0, return_inverse=True)
    W = np.bincount(inv.ravel(), weights=Wall, minlength=len(E))
    med = np.median(np.abs(W)) + 1e-12
    return E, np.clip(W, 0.02 * med, 50.0 * med)


def minimal_surface(V, F, fixed, outer=4, cg_iter=400):
    """Pinkall-Polthier iteration: repeatedly solve the cotangent-Laplace equation
    with the knot fixed -> converges to a discrete minimal surface (soap film)."""
    V = V.copy()
    nv = len(V)
    free = ~fixed
    for _ in range(outer):
        E, W = cotan_edges(V, F)
        dw = np.bincount(E[:, 0], weights=W, minlength=nv) + \
             np.bincount(E[:, 1], weights=W, minlength=nv)

        def L(x):
            s = np.bincount(E[:, 0], weights=W * x[E[:, 1]], minlength=nv) + \
                np.bincount(E[:, 1], weights=W * x[E[:, 0]], minlength=nv)
            return dw * x - s

        for c in range(3):
            xfix = np.where(free, 0.0, V[:, c])
            b = np.where(free, -L(xfix), 0.0)

            def A(y):
                return np.where(free, L(np.where(free, y, 0.0)), 0.0)

            x = np.where(free, V[:, c], 0.0)          # warm start
            r = b - A(x)
            p = r.copy()
            rs = r @ r
            for _ in range(cg_iter):
                if rs < 1e-24:
                    break
                Ap = A(p)
                al = rs / max(p @ Ap, 1e-30)
                x += al * p
                r -= al * Ap
                rn = r @ r
                p = r + (rn / rs) * p
                rs = rn
            V[:, c] = np.where(free, x, V[:, c])
    return V


def _edge_faces(F):
    """dict (a,b) sorted -> list of face indices"""
    from collections import defaultdict
    em = defaultdict(list)
    for fi, f in enumerate(F):
        for k in range(3):
            a, b = int(f[k]), int(f[(k + 1) % 3])
            em[(a, b) if a < b else (b, a)].append(fi)
    return em


def _has_directed(f, a, b):
    return (f[0] == a and f[1] == b) or (f[1] == a and f[2] == b) or (f[2] == a and f[0] == b)


def orient_structured(F, labels):
    """Orient every membrane and every band as a whole, then glue them along the
    Tait graph (membrane -> band -> membrane). The unavoidable orientation seam of a
    one-sided surface therefore lies exactly on a band/membrane arc instead of
    zig-zagging through the triangles. Returns (faces, orientable?)."""
    from collections import deque
    F = F.copy()
    nf = len(F)
    em = _edge_faces(F)

    # --- components: bands by label, membranes by connectivity -----------------
    comp = np.full(nf, -1, np.int64)
    nb = int(labels.max()) + 1 if len(labels) and labels.max() >= 0 else 0
    band = labels >= 0
    comp[band] = labels[band]
    nxt = nb
    for s0 in range(nf):
        if comp[s0] >= 0:
            continue
        comp[s0] = nxt
        dq = deque([s0])
        while dq:
            fi = dq.popleft()
            f = F[fi]
            for k in range(3):
                a, b = int(f[k]), int(f[(k + 1) % 3])
                for g in em[(a, b) if a < b else (b, a)]:
                    if comp[g] < 0 and not band[g]:
                        comp[g] = nxt
                        dq.append(g)
        nxt += 1
    ncomp = nxt

    # --- orient each component internally ----------------------------------------
    seen = np.zeros(nf, bool)
    for s0 in range(nf):
        if seen[s0]:
            continue
        seen[s0] = True
        dq = deque([s0])
        while dq:
            fi = dq.popleft()
            f = F[fi]
            for k in range(3):
                a, b = int(f[k]), int(f[(k + 1) % 3])
                for g in em[(a, b) if a < b else (b, a)]:
                    if g == fi or seen[g] or comp[g] != comp[fi]:
                        continue
                    if not _has_directed(F[g], b, a):
                        F[g] = F[g][::-1].copy()
                    seen[g] = True
                    dq.append(g)

    # --- component graph with "agree" flags ---------------------------------------
    from collections import defaultdict
    cadj = defaultdict(dict)
    for (a, b), fs in em.items():
        if len(fs) != 2:
            continue
        f, g = fs
        cf, cg = comp[f], comp[g]
        if cf == cg:
            continue
        fa = F[f]
        if _has_directed(fa, a, b):
            agree = _has_directed(F[g], b, a)
        else:
            agree = _has_directed(F[g], a, b)
        if cg not in cadj[cf]:
            cadj[cf][cg] = agree
            cadj[cg][cf] = agree

    # --- BFS over components, biggest membrane first -------------------------------
    flip = np.zeros(ncomp, bool)
    done = np.zeros(ncomp, bool)
    sizes = np.bincount(comp, minlength=ncomp)
    order = [c for c in np.argsort(-sizes) if c >= nb] + list(range(nb))
    for root in order:
        if done[root]:
            continue
        done[root] = True
        dq = deque([root])
        while dq:
            c = dq.popleft()
            for d, agree in cadj[c].items():
                if done[d]:
                    continue
                # d must be flipped if (agree XOR flip[c]) is False
                flip[d] = (not agree) != flip[c]
                done[d] = True
                dq.append(d)
    fl = flip[comp]
    F[fl] = F[fl][:, ::-1]

    # --- orientability check -----------------------------------------------------------
    orientable = True
    for (a, b), fs in em.items():
        if len(fs) == 2:
            f, g = F[fs[0]], F[fs[1]]
            if _has_directed(f, a, b) == _has_directed(g, a, b):
                orientable = False
                break
    return F, orientable


def corner_signs(F, nv):
    """sigma[f, k] = +1/-1 so that sigma * vertex_normal matches face f's own normal.
    Only vertices touching the orientation seam get non-trivial signs."""
    from collections import deque, defaultdict
    em = _edge_faces(F)
    sigma = np.ones(F.shape, np.int8)
    seam = set()
    for (a, b), fs in em.items():
        if len(fs) == 2 and _has_directed(F[fs[0]], a, b) == _has_directed(F[fs[1]], a, b):
            seam.add((a, b))
    if not seam:
        return sigma
    sv = set()
    for a, b in seam:
        sv.add(a); sv.add(b)
    vf = defaultdict(list)
    for fi, f in enumerate(F):
        for k in range(3):
            if int(f[k]) in sv:
                vf[int(f[k])].append((fi, k))
    for v in sv:
        fan = vf[v]
        kof = {fi: k for fi, k in fan}
        sg = {fan[0][0]: 1}
        dq = deque([fan[0][0]])
        while dq:
            fi = dq.popleft()
            f = F[fi]
            for w in (int(x) for x in f if int(x) != v):
                key = (v, w) if v < w else (w, v)
                for g in em[key]:
                    if g == fi or g in sg or g not in kof:
                        continue
                    sg[g] = sg[fi] * (-1 if key in seam else 1)
                    dq.append(g)
        for fi, k in fan:
            sigma[fi, k] = sg.get(fi, 1)
    return sigma


def vertex_normals(V, F, sigma):
    tri = V[F]
    fn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])   # area-weighted
    N = np.zeros_like(V)
    for k in range(3):
        np.add.at(N, F[:, k], fn * sigma[:, k:k + 1])
    N /= np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-12)
    return N


def thicken(V, F, sigma, N, half):
    """Closed thick shell = orientation double cover of the (possibly one-sided)
    surface. Top copy of face f sits at +half along f's own normal, bottom copy at
    -half; the two sheets are stitched along the knot. For a one-sided surface both
    sheets join across the seam into ONE seamless orientable shell."""
    nv = len(V)
    VV = np.vstack([V + half * N, V - half * N])
    cid = np.where(sigma > 0, F, F + nv)          # corner -> top copy
    cbot = np.where(sigma > 0, F + nv, F)         # corner -> bottom copy
    faces = [cid, cbot[:, ::-1]]
    em = _edge_faces(F)
    rim = []
    for (a, b), fs in em.items():
        if len(fs) != 1:
            continue
        fi = fs[0]
        f = F[fi]
        ka = int(np.where(f == a)[0][0])
        kb = int(np.where(f == b)[0][0])
        if not _has_directed(f, a, b):
            ka, kb = kb, ka
        ta, tb = cid[fi, ka], cid[fi, kb]
        ba, bb = cbot[fi, ka], cbot[fi, kb]
        rim.append((tb, ta, ba, bb))
    return VV, [tuple(int(i) for i in t) for t in np.vstack(faces)] + \
               [tuple(int(i) for i in q) for q in rim]


def orient_faces(F):
    """BFS orientation propagation. Returns (faces, orientable?)."""
    F = F.copy()
    from collections import defaultdict, deque
    emap = defaultdict(list)
    for fi, f in enumerate(F):
        for k in range(3):
            a, b = f[k], f[(k + 1) % 3]
            emap[(min(a, b), max(a, b))].append(fi)
    visited = np.zeros(len(F), bool)
    orientable = True

    def directed(f, a, b):
        for k in range(3):
            if f[k] == a and f[(k + 1) % 3] == b:
                return True
        return False

    for s in range(len(F)):
        if visited[s]:
            continue
        visited[s] = True
        dq = deque([s])
        while dq:
            fi = dq.popleft()
            f = F[fi]
            for k in range(3):
                a, b = f[k], f[(k + 1) % 3]
                for nj in emap[(min(a, b), max(a, b))]:
                    if nj == fi:
                        continue
                    # neighbour must traverse edge as b->a
                    ok = directed(F[nj], b, a)
                    if not visited[nj]:
                        if not ok:
                            F[nj] = F[nj][::-1].copy()
                        visited[nj] = True
                        dq.append(nj)
                    elif not ok:
                        orientable = False
    return F, orientable


def boundary_loops(F):
    from collections import defaultdict
    cnt = defaultdict(int)
    for f in F:
        for k in range(3):
            a, b = int(f[k]), int(f[(k + 1) % 3])
            cnt[(min(a, b), max(a, b))] += 1
    adj = defaultdict(list)
    for (a, b), c in cnt.items():
        if c == 1:
            adj[a].append(b)
            adj[b].append(a)
    seen = set()
    loops = []
    for s in list(adj.keys()):
        if s in seen:
            continue
        loop = [s]
        seen.add(s)
        prev, cur = None, s
        while True:
            nxt = [v for v in adj[cur] if v != prev and (v not in seen or v == s)]
            if not nxt:
                break
            v = nxt[0]
            if v == s:
                break
            loop.append(v)
            seen.add(v)
            prev, cur = cur, v
        loops.append(loop)
    return loops


# =============================================================================
#  4. PRESETS  (family + parameters)
# =============================================================================

PRESETS = {
    'MOBIUS_TREFOIL': ("Möbius Trefoil", "Trefoil diagram, lobes shaded: a Möbius band with 3 half-twists",
        dict(family='ROSETTE', windings=2, lobes=3, amp=0.55, shade='A', canvas='PLANE')),
    'TRI_WEB': ("Tri-Web", "3-fold web with three holes and twisted corners",
        dict(family='EPI', f1=1, f2=4, f3=0, a2=0.7, a3=0.0, shade='A', canvas='PLANE')),
    'TRI_SWIRL': ("Tri-Web Swirl", "Tri-Web with a third circle: Shape Phase adds a chiral swirl",
        dict(family='EPI', f1=1, f2=4, f3=7, a2=0.7, a3=0.2, ph3=1.6, shade='A', canvas='PLANE')),
    'TRI_WEAVE_BOWL': ("Tri-Weave Bowl", "3-fold interwoven bands on a spherical bowl",
        dict(family='EPI', f1=1, f2=4, f3=0, a2=1.3, a3=0.0, shade='B', canvas='SPHERE', coverage=0.45)),
    'BORROMEAN': ("Borromean Web", "Three alternating rings, one-sided spanning surface",
        dict(family='RINGS', rings=3, ring_dist=0.55, ring_radius=1.0, ring_squash=1.0, shade='A', canvas='PLANE')),
    'HEXA_STAR': ("Hexa Star", "6-fold star web with diamond twists",
        dict(family='EPI', f1=1, f2=7, f3=0, a2=0.7, a3=0.0, shade='A', canvas='PLANE')),
    'HEXA_SHELL': ("Hexa Shell", "6-fold lace wrapped on a sphere",
        dict(family='EPI', f1=1, f2=7, f3=0, a2=0.7, a3=0.0, shade='A', canvas='SPHERE', coverage=0.45)),
    'QUAD_LACE': ("Quad Lace Shell", "4-fold lace with nested holes",
        dict(family='EPI', f1=1, f2=5, f3=0, a2=1.3, a3=0.0, shade='A', canvas='SPHERE', coverage=0.45)),
    'PENTA_WEAVE': ("Penta Weave Ball", "5-lobed Turk's-head weave",
        dict(family='ROSETTE', windings=3, lobes=5, amp=0.45, shade='A', canvas='SPHERE', coverage=0.5)),
    'TETRA_BALL': ("Tetra Ribbon Ball", "Ribbons braided around a ball",
        dict(family='ROSETTE', windings=4, lobes=3, amp=0.6, shade='A', canvas='SPHERE', coverage=0.5)),
    'PIERCED_DISC': ("Pierced Disc", "Dense ring of elliptical holes on a dome",
        dict(family='ROSETTE', windings=2, lobes=11, amp=0.3, shade='B', canvas='SPHERE', coverage=0.4)),
    'LISSA_LATTICE': ("Lissajous Lattice", "3:4 Lissajous weave on a sphere, 11 crosscaps",
        dict(family='LISSAJOUS', lx=3, ly=4, lphase=0.25, laspect=1.0, shade='B', canvas='SPHERE', coverage=0.45)),
    'RING_CHAIN': ("Ring Chain (two-sided)", "Chain of rings - an orientable counter-example",
        dict(family='CHAIN', rings=4, chain_spacing=1.35, ring_squash=0.7, shade='A', canvas='PLANE')),
    'CUSTOM': ("Custom", "Keep the parameters below as they are", {}),
}


_SUPPRESS = [False]


def apply_preset(self, context):
    key = self.preset
    if key == 'CUSTOM':
        return
    _SUPPRESS[0] = True
    try:
        for k, v in PRESETS[key][2].items():
            setattr(self, k, v)
    finally:
        _SUPPRESS[0] = False


# =============================================================================
#  5. BLENDER SIDE
# =============================================================================

def props_to_dict(pr):
    keys = ["family", "windings", "lobes", "amp", "f1", "f2", "f3", "a2", "a3", "ph2",
            "ph3", "lx", "ly", "lphase", "laspect", "rings", "ring_dist", "ring_radius",
            "ring_squash", "chain_spacing", "shade", "crossing_rule", "handedness", "seed",
            "band_size", "twist_height", "band_res", "density", "inflate", "relax_iter",
            "soap_film", "soap_iter", "rotation",
            "dome", "size", "strand_flatten", "canvas", "coverage", "edge_smooth"]
    return {k: getattr(pr, k) for k in keys}


def get_material(name, base, dark):
    mat = bpy.data.materials.get(name)
    if mat:
        return mat
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial"); out.location = (500, 0)
    bs = nt.nodes.new("ShaderNodeBsdfPrincipled"); bs.location = (200, 0)
    bs.inputs["Roughness"].default_value = 0.45
    lw = nt.nodes.new("ShaderNodeLayerWeight"); lw.location = (-400, 100)
    lw.inputs["Blend"].default_value = 0.35
    cr = nt.nodes.new("ShaderNodeValToRGB"); cr.location = (-150, 100)
    cr.color_ramp.elements[0].color = base
    cr.color_ramp.elements[1].color = dark
    nt.links.new(lw.outputs["Facing"], cr.inputs["Fac"])
    nt.links.new(cr.outputs["Color"], bs.inputs["Base Color"])
    nt.links.new(bs.outputs["BSDF"], out.inputs["Surface"])
    return mat


def _ensure_outward(me):
    """Flip the whole closed shell if its signed volume is negative."""
    vol = 0.0
    co = [v.co for v in me.vertices]
    for p in me.polygons:
        vs = p.vertices
        a = co[vs[0]]
        for i in range(1, len(vs) - 1):
            vol += a.dot(co[vs[i]].cross(co[vs[i + 1]]))
    if vol < 0:
        me.flip_normals()


def make_objects(context, V, F, loops, pr):
    for nm in (OBJ_NAME, EDGE_NAME):
        o = bpy.data.objects.get(nm)
        if o:
            data = o.data
            bpy.data.objects.remove(o, do_unlink=True)
            if data and data.users == 0:
                if isinstance(data, bpy.types.Mesh):
                    bpy.data.meshes.remove(data)
                else:
                    bpy.data.curves.remove(data)

    F = np.asarray(F, dtype=np.int64)
    sigma = corner_signs(F, len(V))
    N = vertex_normals(V, F, sigma)
    me = bpy.data.meshes.new(OBJ_NAME)
    if pr.use_solidify and pr.thickness > 0:
        # seamless closed shell (orientation double cover) instead of the modifier
        half = 0.5 * pr.thickness * pr.size * 0.01
        VV, FF = thicken(V, F, sigma, N, half)
        me.from_pydata([tuple(v) for v in VV], [], FF)
        me.validate()
        me.shade_smooth() if pr.smooth else me.shade_flat()
        _ensure_outward(me)
    else:
        me.from_pydata([tuple(v) for v in V], [], [tuple(int(i) for i in f) for f in F])
        ok = not me.validate()
        if pr.smooth:
            me.shade_smooth()
            if ok and len(me.polygons) == len(F):
                # per-corner normals that always agree with their own face:
                # no dark seam on the one-sided surface
                ln = (N[F] * sigma[:, :, None]).reshape(-1, 3)
                me.normals_split_custom_set([tuple(n) for n in ln])
        else:
            me.shade_flat()
    obj = bpy.data.objects.new(OBJ_NAME, me)
    context.collection.objects.link(obj)
    if pr.subdiv > 0:
        m = obj.modifiers.new("Subdivision", 'SUBSURF')
        m.levels = pr.subdiv
        m.render_levels = pr.subdiv
        if hasattr(m, "use_custom_normals"):
            m.use_custom_normals = True
    if pr.add_material:
        me.materials.append(get_material(MAT_NAME, (0.80, 0.90, 0.92, 1), (0.42, 0.62, 0.68, 1)))

    if pr.edge_tube > 0:
        cu = bpy.data.curves.new(EDGE_NAME, 'CURVE')
        cu.dimensions = '3D'
        cu.bevel_depth = pr.edge_tube * pr.size * 0.004
        cu.bevel_resolution = 2
        for lp in loops:
            if len(lp) < 3:
                continue
            sp = cu.splines.new('POLY')
            sp.points.add(len(lp) - 1)
            for k, vi in enumerate(lp):
                x, y, z = V[vi]
                sp.points[k].co = (x, y, z, 1)
            sp.use_cyclic_u = True
        eo = bpy.data.objects.new(EDGE_NAME, cu)
        context.collection.objects.link(eo)
        eo.parent = obj
        if pr.add_material:
            cu.materials.append(get_material(EDGE_MAT_NAME, (0.25, 0.35, 0.38, 1), (0.1, 0.15, 0.17, 1)))

    for o in context.selected_objects:
        o.select_set(False)
    obj.select_set(True)
    context.view_layer.objects.active = obj
    return obj


def live_update(self, context):
    if self.live and not _SUPPRESS[0]:
        bpy.ops.mesh.metamobius_generate()


FAMILIES = [
    ('ROSETTE', "Rosette / Turk's Head", "r = 1 + A·cos(L t), angle = W t"),
    ('LISSAJOUS', "Lissajous", "x = sin(a t + φ), y = sin(b t)"),
    ('EPI', "Epitrochoid / Farris Wheel", "z = e^(i f1 t) + a2 e^(i f2 t) + a3 e^(i f3 t) - (f2-f1)-fold symmetric"),
    ('RINGS', "Ring Link (Borromean …)", "m rings arranged in a circle"),
    ('CHAIN', "Ring Chain", "m rings in a row"),
]


def preset_update(self, context):
    apply_preset(self, context)
    live_update(self, context)


class MetamobiusProps(PropertyGroup):
    preset: EnumProperty(name="Surface", items=[(k, v[0], v[1]) for k, v in PRESETS.items()],
                         default='MOBIUS_TREFOIL', update=preset_update)
    live: BoolProperty(name="Live Update", description="Regenerate automatically whenever a parameter changes",
                       default=False, update=lambda s, c: live_update(s, c))

    family: EnumProperty(name="Diagram", items=FAMILIES, default='ROSETTE', update=live_update)
    windings: IntProperty(name="Windings W", default=2, min=1, max=15, update=live_update)
    lobes: IntProperty(name="Lobes L", default=3, min=1, max=31, update=live_update)
    amp: FloatProperty(name="Amplitude", default=0.55, min=0.02, max=2.5, update=live_update)

    f1: IntProperty(name="Freq 1", default=1, min=-20, max=20, update=live_update)
    f2: IntProperty(name="Freq 2", default=-5, min=-20, max=20, update=live_update)
    f3: IntProperty(name="Freq 3", default=7, min=-20, max=20, update=live_update)
    a2: FloatProperty(name="Amp 2", default=0.55, min=0.0, max=2.0, update=live_update)
    a3: FloatProperty(name="Amp 3", default=0.25, min=0.0, max=2.0, update=live_update)
    ph2: FloatProperty(name="Phase 2", default=0.0, subtype='ANGLE', update=live_update)
    ph3: FloatProperty(name="Shape Phase", description="Phase of the 3rd circle - only acts when Amp 3 > 0 "
                       "(breaks mirror symmetry, adds a chiral swirl)",
                       default=0.5, subtype='ANGLE', update=live_update)
    rotation: FloatProperty(name="Rotation", description="Rotate the diagram (does not change the shape)",
                            default=0.0, subtype='ANGLE', update=live_update)

    lx: IntProperty(name="Freq X", default=3, min=1, max=15, update=live_update)
    ly: IntProperty(name="Freq Y", default=2, min=1, max=15, update=live_update)
    lphase: FloatProperty(name="Phase", default=0.35, subtype='ANGLE', update=live_update)
    laspect: FloatProperty(name="Aspect Y", default=0.75, min=0.2, max=2.0, update=live_update)

    rings: IntProperty(name="Rings", default=3, min=2, max=12, update=live_update)
    ring_dist: FloatProperty(name="Center Distance", default=0.55, min=0.05, max=3.0, update=live_update)
    ring_radius: FloatProperty(name="Ring Radius", default=1.0, min=0.1, max=3.0, update=live_update)
    ring_squash: FloatProperty(name="Ring Squash", default=1.0, min=0.2, max=1.5, update=live_update)
    chain_spacing: FloatProperty(name="Spacing", default=1.35, min=0.3, max=1.95, update=live_update)

    shade: EnumProperty(name="Colouring", items=[
        ('A', "Colour A", "Bounded regions with odd winding become membranes"),
        ('B', "Colour B", "The complementary regions (diagram inverted so they are bounded)")],
        default='A', update=live_update)
    crossing_rule: EnumProperty(name="Crossings", items=[
        ('ALTERNATING', "Alternating (Weave)", "Over/under alternates - all bands twist the same way"),
        ('TORUS', "Torus Knot (Rosette only)", "Crossings from the real (W,L) torus knot"),
        ('RANDOM', "Random", "Random over/under (seeded)")],
        default='ALTERNATING', update=live_update)
    handedness: BoolProperty(name="Mirror Twist", default=False, update=live_update)
    seed: IntProperty(name="Seed", default=1, update=live_update)

    band_size: FloatProperty(name="Band Size", default=1.3, min=0.2, max=3.0, update=live_update)
    twist_height: FloatProperty(name="Twist Height", default=1.5, min=0.0, max=5.0, update=live_update)
    band_res: IntProperty(name="Band Resolution", default=8, min=4, max=32, update=live_update)
    density: IntProperty(name="Mesh Density", default=70, min=20, max=250, update=live_update)
    inflate: FloatProperty(name="Inflate Membranes", default=0.0, min=-5.0, max=5.0, update=live_update)
    relax_iter: IntProperty(name="(legacy)", default=40, min=0, max=500)
    soap_film: FloatProperty(name="Soap Film", description="Blend towards the minimal surface spanning the knot "
                             "(0 = original, 1 = true soap film)", default=0.0, min=0.0, max=1.0,
                             subtype='FACTOR', update=live_update)
    soap_iter: IntProperty(name="Soap Iterations", description="Minimal-surface solver passes "
                           "(more = closer to the true soap film, slower)",
                           default=4, min=1, max=20, update=live_update)
    dome: FloatProperty(name="Dome", default=0.0, min=-2.0, max=2.0, update=live_update)
    strand_flatten: FloatProperty(name="Strand Flatten", description="Pull the knot towards the canvas between crossings",
                                  default=0.3, min=0.0, max=1.0, update=live_update)
    edge_smooth: IntProperty(name="Edge Smoothing", default=6, min=0, max=100, update=live_update)
    canvas: EnumProperty(name="Canvas", items=[
        ('PLANE', "Plane", "Diagram drawn on a plane"),
        ('SPHERE', "Sphere", "Diagram wrapped onto a sphere (stereographic) - shells, bowls, balls")],
        default='PLANE', update=live_update)
    coverage: FloatProperty(name="Sphere Coverage", description="How far the diagram reaches around the sphere",
                            default=0.45, min=0.1, max=0.95, update=live_update)
    size: FloatProperty(name="Size", default=2.0, min=0.01, max=100.0, update=live_update)

    use_solidify: BoolProperty(name="Solidify", description="Closed, seamless thick shell "
                               "(built directly, not with the Solidify modifier)",
                               default=True, update=live_update)
    thickness: FloatProperty(name="Thickness", default=0.6, min=0.01, max=10.0, update=live_update)
    subdiv: IntProperty(name="Subdivision", default=0, min=0, max=3, update=live_update)
    smooth: BoolProperty(name="Shade Smooth", default=True, update=live_update)
    edge_tube: FloatProperty(name="Edge Outline", default=1.0, min=0.0, max=10.0, update=live_update)
    add_material: BoolProperty(name="Preview Material", default=True, update=live_update)

    info: StringProperty(default="")


class METAMOBIUS_OT_generate(Operator):
    bl_idname = "mesh.metamobius_generate"
    bl_label = "Generate Surface"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        pr = context.scene.metamobius
        try:
            V, F, loops, info = build_surface(props_to_dict(pr), report=lambda m: self.report({'INFO'}, m))
        except Exception as ex:
            self.report({'ERROR'}, str(ex))
            pr.info = "Error: " + str(ex)
            return {'CANCELLED'}
        make_objects(context, V, F, loops, pr)
        chi, b = info["chi"], info["boundaries"]
        if info["orientable"]:
            gen = (2 - chi - b) / 2
            topo = f"two-sided | genus {gen:g}"
        else:
            topo = f"ONE-SIDED | {2 - chi - b} crosscaps"
        pr.info = (f"{topo} | χ={chi} | boundary loops={b} | crossings={info['crossings']}")
        return {'FINISHED'}


class METAMOBIUS_PT_panel(Panel):
    bl_label = "Metamobius Surfaces"
    bl_idname = "METAMOBIUS_PT_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Metamobius"

    def draw(self, context):
        pr = context.scene.metamobius
        L = self.layout
        L.prop(pr, "preset")
        row = L.row(align=True)
        row.scale_y = 1.4
        row.operator("mesh.metamobius_generate", icon='MOD_CAST')
        row = L.row(align=True)
        row.prop(pr, "live", text="Live Update  ON" if pr.live else "Live Update  OFF",
                 icon='CHECKBOX_HLT' if pr.live else 'CHECKBOX_DEHLT', toggle=True)
        b = L.box()
        b.label(text="Diagram", icon='CURVE_BEZCIRCLE')
        b.prop(pr, "family", text="")
        f = pr.family
        if f == 'ROSETTE':
            b.prop(pr, "windings"); b.prop(pr, "lobes"); b.prop(pr, "amp")
        elif f == 'EPI':
            c = b.column(align=True)
            c.prop(pr, "f1"); c.prop(pr, "f2"); c.prop(pr, "f3")
            c = b.column(align=True)
            c.prop(pr, "a2"); c.prop(pr, "a3")
            if pr.a3 > 1e-6:
                c.prop(pr, "ph3")
            else:
                c.label(text="Amp 3 > 0 enables Shape Phase", icon='INFO')
        elif f == 'LISSAJOUS':
            b.prop(pr, "lx"); b.prop(pr, "ly"); b.prop(pr, "lphase"); b.prop(pr, "laspect")
        elif f == 'RINGS':
            b.prop(pr, "rings"); b.prop(pr, "ring_dist"); b.prop(pr, "ring_radius"); b.prop(pr, "ring_squash")
        elif f == 'CHAIN':
            b.prop(pr, "rings"); b.prop(pr, "chain_spacing"); b.prop(pr, "ring_squash")

        b.prop(pr, "rotation")

        b = L.box()
        b.label(text="Surface", icon='MOD_WIREFRAME')
        b.prop(pr, "shade")
        r = b.row(align=True)
        r.prop(pr, "canvas", expand=True)
        if pr.canvas == 'SPHERE':
            b.prop(pr, "coverage")
        b.prop(pr, "crossing_rule")
        r = b.row(); r.prop(pr, "handedness")
        if pr.crossing_rule == 'RANDOM':
            r.prop(pr, "seed")
        c = b.column(align=True)
        c.prop(pr, "band_size"); c.prop(pr, "twist_height"); c.prop(pr, "band_res")
        c = b.column(align=True)
        c.prop(pr, "inflate"); c.prop(pr, "edge_smooth")
        c = b.column(align=True)
        c.prop(pr, "soap_film")
        sub = c.row(align=True)
        sub.active = pr.soap_film > 0
        sub.prop(pr, "soap_iter")
        c.prop(pr, "strand_flatten")
        if pr.canvas == 'PLANE':
            c.prop(pr, "dome")
        b.prop(pr, "density")

        b = L.box()
        b.label(text="Output", icon='SHADING_RENDERED')
        b.prop(pr, "size")
        r = b.row(align=True)
        r.prop(pr, "use_solidify")
        sub = r.row(align=True)
        sub.active = pr.use_solidify
        sub.prop(pr, "thickness")
        b.prop(pr, "subdiv")
        b.prop(pr, "edge_tube"); b.prop(pr, "smooth"); b.prop(pr, "add_material")


classes = (MetamobiusProps, METAMOBIUS_OT_generate, METAMOBIUS_PT_panel)


def register():
    for c in classes:
        bpy.utils.register_class(c)
    bpy.types.Scene.metamobius = PointerProperty(type=MetamobiusProps)


def unregister():
    del bpy.types.Scene.metamobius
    for c in reversed(classes):
        bpy.utils.unregister_class(c)


if __name__ == "__main__":
    register()
