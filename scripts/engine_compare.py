#!/usr/bin/env python3
"""Both engines on the same models, as SVG plots and one markdown page.

Writes docs/results/engine_comparison.md and docs/results/plots/*.svg.

For every kernel model in examples/ (the GSSK schema):
  - the kernel is run under RK4, Euler and its expm method, at the model's dt;
  - where the projection carries the model whole, the reference is the
    Giannantoni engine's exp(A t) at the same times (bin/gia_bridge --csv);
  - otherwise there is no exact solution to compare with, and the reference is
    the kernel's own RK4 at dt/8, labelled as such, beside the coverage
    report's reasons the engine cannot carry the model.
  The plot shows each storage's reference trajectory with the model's own
  method on top, and the largest scaled difference over time for each method.

For every Giannantoni-format model in examples/giannantoni/, which the kernel
cannot read: the engine's network trajectories and its output-projection
drift, [09 Eq 13], the one incipient quantity it reports per component.

Both engines' network trajectories are classical (PLAN.md §2 E3), so every
difference in the kernel plots is integration error, never incipient drift.

    python3 -I scripts/engine_compare.py            # rewrite the page and plots
    python3 -I scripts/engine_compare.py --check    # fail if they would change

Verifies: FR-OUT-005 (T-OUT-05)
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
OUT_MD = os.path.join(ROOT, "docs", "results", "engine_comparison.md")
OUT_DIR = os.path.join(ROOT, "docs", "results", "plots")
METHODS = ["rk4", "euler", "expm"]
FLOOR = 1e-13          # differences below this are rounding; drawn on the floor
FINE = 8               # reference refinement when there is no exact solution
FINER = 64             # ... and the run that says how good that reference is
MAX_PANELS = 4         # storages drawn as small multiples; the rest are tabled
GIA_STEPS = 100
MAX_POINTS = 400       # polyline vertices per series

# Reference palette (dataviz skill, references/palette.md), validated with
# scripts/validate_palette.js: slots 1-3 all-pairs and slots 1-8 adjacent, in
# both modes. Light slots 3-5 sit below 3:1, so every series is direct-labelled
# and every plot has a table of the same numbers beside it.
SLOTS_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SLOTS_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]
METHOD_SLOT = {m: i for i, m in enumerate(METHODS)}   # colour follows the method, everywhere


# ---------------------------------------------------------------- data

def read_csv(path):
    with open(path, newline="") as f:
        rows = [r for r in csv.reader(f) if r]
    head = rows[0]
    cols = {h: [] for h in head}
    for r in rows[1:]:
        for h, v in zip(head, r):
            cols[h].append(float(v))
    return head, cols


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def kernel_run(doc, tmp, tag, method, dt=None):
    d = json.loads(json.dumps(doc))
    d.setdefault("config", {})["method"] = method
    if dt is not None:
        d["config"]["dt"] = dt
    mp, cp = os.path.join(tmp, tag + ".json"), os.path.join(tmp, tag + ".csv")
    with open(mp, "w") as f:
        json.dump(d, f)
    r = run([os.path.join(BIN, "gssk"), mp, cp, "--full-precision"])
    if r.returncode != 0 or not os.path.exists(cp):
        return None
    return read_csv(cp)[1]


def bridge_run(model, method, tmp):
    cp = os.path.join(tmp, "bridge_%s.csv" % method)
    r = run([os.path.join(BIN, "gia_bridge"), model, "--method", method, "--csv", cp])
    cov, reasons = None, []
    for line in (r.stdout + r.stderr).splitlines():
        s = line.strip()
        if s.startswith("bridge.coverage:"):
            cov = s.split(":", 1)[1].strip()
        elif s.startswith("- "):
            reasons.append(s[2:])
    if r.returncode != 0:
        return None, cov, reasons
    return read_csv(cp)[1], cov, reasons


def scaled_error(ref, got, nodes):
    """Per time step, the largest |got - ref| over nodes, each scaled by that
    node's largest |ref| over the run (the survey's metric, FR-BRG-002)."""
    n = min(len(ref[nodes[0]]), len(got[nodes[0]]))
    out = [0.0] * n
    for k in nodes:
        mag = max(abs(v) for v in ref[k]) or 1.0
        for i in range(n):
            out[i] = max(out[i], abs(got[k][i] - ref[k][i]) / mag)
    return out


# ---------------------------------------------------------------- svg

def fmt(v):
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e4 or a < 1e-2:
        m, e = ("%.0e" % v).split("e")
        return "%se%d" % (m, int(e))
    return ("%.3g" % v)


def nice_ticks(lo, hi, n=4):
    if hi <= lo:
        hi = lo + (abs(lo) or 1.0)
    span = hi - lo
    raw = span / n
    mag = 10 ** __import__("math").floor(__import__("math").log10(raw))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), default=raw)
    import math
    t0 = math.ceil(lo / step) * step
    ticks, t = [], t0
    while t <= hi + step * 1e-9:
        ticks.append(0.0 if abs(t) < step * 1e-9 else t)
        t += step
    return ticks


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


STYLE = """<style>
  .bg { fill: #fcfcfb; } .ink { fill: #0b0b0b; } .ink2 { fill: #52514e; } .muted { fill: #898781; }
  .grid { stroke: #e1e0d9; stroke-width: 1; } .axis { stroke: #c3c2b7; stroke-width: 1; }
  .ref { stroke: #0b0b0b; stroke-width: 1.5; fill: none; }
  text { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; }
%s
  @media (prefers-color-scheme: dark) {
    .bg { fill: #1a1a19; } .ink { fill: #f0efec; } .ink2 { fill: #c3c2b7; }
    .grid { stroke: #2c2c2a; } .axis { stroke: #383835; } .ref { stroke: #f0efec; }
%s
  }
</style>"""


def slot_css():
    light = "\n".join("  .s%d { stroke: %s; } .f%d { fill: %s; }" % (i, c, i, c) for i, c in enumerate(SLOTS_LIGHT))
    dark = "\n".join("    .s%d { stroke: %s; } .f%d { fill: %s; }" % (i, c, i, c) for i, c in enumerate(SLOTS_DARK))
    return STYLE % (light, dark)


class Svg:
    def __init__(self, w, h, title, desc):
        self.w, self.h, self.parts = w, h, []
        self.head = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
                     'role="img" aria-labelledby="t d">\n<title id="t">%s</title>\n<desc id="d">%s</desc>\n%s\n'
                     '<rect class="bg" x="0" y="0" width="%d" height="%d" rx="8"/>\n'
                     % (w, h, w, h, esc(title), esc(desc), slot_css(), w, h))

    def add(self, s):
        self.parts.append(s)

    def text(self, x, y, s, cls="ink2", size=12, anchor="start", weight=None):
        wt = ' font-weight="%s"' % weight if weight else ""
        self.add('<text class="%s" x="%.1f" y="%.1f" font-size="%d" text-anchor="%s"%s>%s</text>'
                 % (cls, x, y, size, anchor, wt, esc(s)))

    def done(self):
        return self.head + "\n".join(self.parts) + "\n</svg>\n"


def panel(svg, x, y, w, h, xs, series, title, logy=False, xlabel="t", right=92):
    """series: list of (label, ys, kind, slot) with kind 'ref' (ink line),
    'line' (slot line) or 'dots' (slot markers). Draws axes, grid, marks and a
    direct label at each line's right end."""
    import math
    L, R, T, B = 52, right, 22, 30
    pw, ph = w - L - R, h - T - B
    svg.text(x + L, y + 14, title, cls="ink", size=12, weight="600")
    allx = xs
    x0, x1 = min(allx), max(allx)
    if x1 == x0:
        x1 = x0 + 1
    vals = [v for _, ys, _, _ in series for v in ys]
    if logy:
        vals = [max(v, FLOOR) for v in vals]
        lo = math.floor(math.log10(min(vals)))
        hi = math.ceil(math.log10(max(vals)))
        if hi <= lo:
            hi = lo + 1
        lo = max(lo, int(math.log10(FLOOR)))
        tf = lambda v: math.log10(max(v, FLOOR))
        ticks = list(range(lo, hi + 1))
        step = max(1, len(ticks) // 5)
        ticks = ticks[::step]
        ylo, yhi = lo, hi
        tick_label = lambda t: "1e%d" % t
        tick_pos = lambda t: t
    else:
        ylo, yhi = min(vals), max(vals)
        if yhi - ylo < 1e-12 * max(1.0, abs(yhi)):
            ylo, yhi = ylo - 1, yhi + 1
        ticks = nice_ticks(ylo, yhi)
        ylo, yhi = min(ylo, ticks[0]), max(yhi, ticks[-1])
        tf = lambda v: v
        tick_label = fmt
        tick_pos = lambda t: t
    px = lambda v: x + L + (v - x0) / (x1 - x0) * pw
    py = lambda v: y + T + ph - (tf(v) - ylo) / (yhi - ylo) * ph
    pyt = lambda t: y + T + ph - (tick_pos(t) - ylo) / (yhi - ylo) * ph
    for t in ticks:
        yy = pyt(t)
        svg.add('<line class="grid" x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f"/>' % (x + L, x + L + pw, yy, yy))
        svg.text(x + L - 6, yy + 4, tick_label(t), cls="muted", size=10, anchor="end")
    svg.add('<line class="axis" x1="%.1f" x2="%.1f" y1="%.1f" y2="%.1f"/>' % (x + L, x + L + pw, y + T + ph, y + T + ph))
    for t in nice_ticks(x0, x1, 4):
        if x0 - 1e-9 <= t <= x1 + 1e-9:
            svg.text(px(t), y + T + ph + 14, fmt(t), cls="muted", size=10, anchor="middle")
    svg.text(x + L + pw, y + T + ph + 26, xlabel, cls="muted", size=10, anchor="end")
    labels = []
    for label, ys, kind, slot in series:
        n = min(len(xs), len(ys))
        if kind == "dots":
            every = max(1, n // 16)
            for i in range(0, n, every):
                svg.add('<circle class="f%d" cx="%.1f" cy="%.1f" r="4" stroke="none"><title>%s t=%s: %s</title></circle>'
                        % (slot, px(xs[i]), py(ys[i]), esc(label), fmt(xs[i]), fmt(ys[i])))
            continue
        if n > MAX_POINTS:
            # One vertex per bucket. On a log (error) panel each keeps its
            # bucket's largest value, so a spike is never thinned away.
            idx, pick = [], []
            for b0 in range(MAX_POINTS):
                lo_i, hi_i = b0 * n // MAX_POINTS, max(b0 * n // MAX_POINTS + 1, (b0 + 1) * n // MAX_POINTS)
                j = max(range(lo_i, hi_i), key=lambda q: ys[q]) if logy else lo_i
                idx.append(j)
            if idx[-1] != n - 1:
                idx.append(n - 1)
        else:
            idx = range(n)
        pts = " ".join("%.1f,%.1f" % (px(xs[i]), py(ys[i])) for i in idx)
        cls = "ref" if kind == "ref" else "s%d" % slot
        extra = "" if kind == "ref" else ' stroke-width="2" fill="none" stroke-linejoin="round"'
        svg.add('<polyline class="%s" points="%s"%s><title>%s</title></polyline>' % (cls, pts, extra, esc(label)))
        labels.append([py(ys[n - 1]), label, "ink" if kind == "ref" else "ink2"])
    # direct labels at the right edge, nudged apart so they never overlap
    labels.sort()
    for i in range(1, len(labels)):
        if labels[i][0] - labels[i - 1][0] < 12:
            labels[i][0] = labels[i - 1][0] + 12
    for yy, label, cls in labels:
        svg.text(x + L + pw + 6, yy + 4, label, cls=cls, size=10)


# ---------------------------------------------------------------- kernel models

def kernel_model(model, tmp):
    name = os.path.splitext(os.path.basename(model))[0]
    with open(model) as f:
        doc = json.load(f)
    cfg = doc.get("config", {})
    own = cfg.get("method", "rk4")
    types = {n.get("id"): n.get("type") for n in doc.get("nodes", [])}

    runs, cov, reasons, exact = {}, None, [], None
    for m in METHODS:
        got, c, why = bridge_run(model, m, tmp)
        cov = cov or c
        reasons = reasons or why
        if got is None:
            break
        ids = [h[:-len("_kernel")] for h in got if h.endswith("_kernel")]
        runs[m] = {k: got[k + "_kernel"] for k in ids}
        if exact is None:
            exact = {k: got[k + "_gia"] for k in ids}
        exact_t = got["time"]
    whole = len(runs) == len(METHODS)

    finer = None
    if whole:
        ref, t, ref_label = exact, exact_t, "Giannantoni exp(At)"
    else:
        runs = {}
        for m in METHODS:
            r = kernel_run(doc, tmp, name + "_" + m, m)
            if r is not None:
                runs[m] = r
        dt = cfg.get("dt", 0.1)
        fine = kernel_run(doc, tmp, name + "_fine", "rk4", dt / FINE)
        if fine is None or "rk4" not in runs:
            return None
        t = runs["rk4"]["time"]
        ref = {k: v[::FINE][:len(t)] for k, v in fine.items()}
        ref_label = "kernel RK4 at dt/%d" % FINE
        finer = kernel_run(doc, tmp, name + "_finer", "rk4", dt / FINER)
    nodes = [k for k in ref if k != "time" and all(k in r for r in runs.values())]
    stores = [k for k in nodes if types.get(k) == "storage"] or nodes
    err = {m: scaled_error(ref, runs[m], stores) for m in runs}
    # How far the reference itself is from converged: exact for exp(A t);
    # otherwise its own difference from a run eight times finer still.
    ref_err = 0.0
    if not whole and finer is not None:
        ref_err = max(scaled_error(ref, {k: v[::FINER][:len(t)] for k, v in finer.items()}, stores) or [0.0])
    worst = {m: max(e) for m, e in err.items()}

    # ---- figure: small multiples of storages, then the error over time
    shown = stores[:MAX_PANELS]
    cols = max(1, len(shown))
    W = 760
    pw = W // cols if cols > 1 else W
    H = 40 + 190 + 230
    svg = Svg(W, H, "%s: kernel against %s" % (name, ref_label),
              "Storages under %s (markers) over the %s reference (line), and the largest scaled "
              "difference from that reference over time for kernel RK4, Euler and expm." % (own, ref_label))
    svg.text(16, 24, name, cls="ink", size=14, weight="700")
    svg.text(16 + 8 * len(name) + 18, 24, "line: %s   markers: kernel %s (the model's own method)"
             % (ref_label, own), cls="ink2", size=11)
    own_slot = METHOD_SLOT.get(own, 0)
    for j, k in enumerate(shown):
        short = "exp(At)" if whole else "RK4 dt/%d" % FINE
        panel(svg, j * pw, 36, pw, 190, t, [(short, ref[k], "ref", 0),
                                            ("kernel " + own, runs.get(own, runs["rk4"])[k], "dots", own_slot)],
              k, right=(60 if cols > 1 else 92))
    ser = [("RK4", [max(v, FLOOR) for v in err["rk4"]], "line", METHOD_SLOT["rk4"])]
    for m, lab in (("euler", "Euler"), ("expm", "expm")):
        if m in err:
            ser.append((lab, [max(v, FLOOR) for v in err[m]], "line", METHOD_SLOT[m]))
    panel(svg, 0, 36 + 196, W, 220, t, ser,
          "Largest scaled |kernel − %s| over the storages  (log scale)" % ref_label, logy=True)
    with open(os.path.join(tmp, "plots", name + ".svg"), "w") as f:
        f.write(svg.done())

    return {"name": name, "own": own, "whole": whole, "cov": cov or "—", "reasons": reasons,
            "ref": ref_label, "ref_err": ref_err, "worst": worst, "stores": stores, "shown": shown,
            "desc": doc.get("metadata", {}).get("description", "")}


# ---------------------------------------------------------------- Giannantoni models

def gia_model(model, tmp):
    name = os.path.splitext(os.path.basename(model))[0]
    cp = os.path.join(tmp, "gia_" + name + ".csv")
    r = run([os.path.join(BIN, "giannantoni_sim"), model, "--csv", cp, "--steps", str(GIA_STEPS),
             "--out", os.path.join(tmp, "gia_" + name + ".out.json")])
    if r.returncode != 0:
        return {"name": name, "error": (r.stderr or r.stdout).strip().splitlines()[-1:]}
    head, cols = read_csv(cp)
    comps = [h[:-2] for h in head if h.endswith("_Q")]
    t = cols["time"]
    with open(model) as f:
        doc = json.load(f)
    title = doc.get("system", doc.get("name", name))
    W, H = 760, 40 + 230 + 230
    svg = Svg(W, H, "%s: Giannantoni engine" % name,
              "Network trajectories Q(t) = exp(At)Q(0) per component, and each component's "
              "output-projection drift [09 Eq 13].")
    svg.text(16, 24, name, cls="ink", size=14, weight="700")
    svg.text(16 + 8 * len(name) + 18, 24, "Giannantoni engine — the kernel cannot read this format",
             cls="ink2", size=11)
    sl = {c: i % len(SLOTS_LIGHT) for i, c in enumerate(comps)}
    panel(svg, 0, 36, W, 226, t, [(c, cols[c + "_Q"], "line", sl[c]) for c in comps],
          "Network trajectories Q(t) = exp(A t) Q(0)  (classical, PLAN E3)")
    panel(svg, 0, 36 + 232, W, 226, t, [(c, cols[c + "_drift_proj"], "line", sl[c]) for c in comps],
          "Output-projection drift per component  (incipient, [09 Eq 13])")
    with open(os.path.join(tmp, "plots", "gia_" + name + ".svg"), "w") as f:
        f.write(svg.done())
    rows = []
    for c in comps:
        q, d = cols[c + "_Q"], cols[c + "_drift_proj"]
        rows.append((c, q[0], q[-1], max(abs(v) for v in d)))
    return {"name": name, "title": title, "rows": rows, "t_end": t[-1]}


# ---------------------------------------------------------------- page

def sci(v):
    return "< %.0e" % FLOOR if v < FLOOR else "%.2e" % v


def page(kern, gia):
    whole = [k for k in kern if k["whole"]]
    out = [
        "# Both engines on the same models",
        "",
        "Generated by `python3 -I scripts/engine_compare.py` (`make engine-comparison`), and checked against "
        "a fresh run in CI (`make check-engine-comparison`, T-OUT-05). Do not edit by hand.",
        "",
        "**What is compared.** Every kernel model in `examples/` is run by the kernel under RK4, Euler and "
        "its `expm` method. Where the projection (ADR 0011) carries the model whole, the reference is the "
        "Giannantoni engine's `exp(A t)` at the same step times, so the difference is each kernel method's "
        "error against the exact solution. Where it does not, there is no exact solution to compare with: "
        "the reference is the kernel's own RK4 at dt/%d, and the coverage report's reasons are listed. "
        "Both engines' network trajectories are classical (PLAN.md §2 E3), so no difference here is "
        "incipient drift. The Giannantoni-format models are shown with the one incipient quantity the "
        "engine reports per component, the output-projection drift." % FINE,
        "",
        "Differences are **scaled**: |Δ| over the node's largest magnitude in the run, maximised over the "
        "storages (the Level 1 survey's metric). Below %.0e is rounding and is drawn on the floor." % FLOOR,
        "",
        "**Reading the plots.** In exponential form RK4 sits near 1e-9 and Euler near 1e-3, and expm is exact "
        "to rounding, because it *is* `exp(A t)` for an autonomous linear model. Three things move these "
        "orders: a nonlinear law, which expm only linearises per step; forcing, which expm holds fixed across "
        "a step; and a jump in a forcing (a `step`, `square` or step-interpolated `table`), which makes every "
        "method first order at the jump, so RK4's advantage disappears there. A fast transient the model's dt "
        "barely resolves (in the Odum countercurrent models, a price that relaxes in a few steps) puts RK4's "
        "largest difference at the start of the run. And a flow through a processing node (interaction, "
        "gain, exchange, switch, loop-limited) enters expm as a constant over the step, evaluated at its "
        "start (`idc_step_ex` in `src/gssk.c`): where such flows drive the model, expm is first order, "
        "and its figure equals Euler's exactly (`atwood_model`, `economic_model`, the `exchange_price_node` "
        "models).",
        "",
        "## Kernel models",
        "",
        "%d of %d are carried whole and compared with `exp(A t)`." % (len(whole), len(kern)),
        "",
        "*Reference error* is how far the reference is from converged: zero for `exp(A t)`, and for "
        "RK4 at dt/%d its own difference from RK4 at dt/%d. A method's figure is only meaningful where it "
        "is well above that." % (FINE, FINER),
        "",
        "| Model | Coverage | Reference | Reference error | Own method | RK4 | Euler | expm |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for k in kern:
        w = k["worst"]
        out.append("| [`%s`](#%s) | %s | %s | %s | %s | %s | %s | %s |" % (
            k["name"], k["name"], k["cov"], k["ref"], "exact" if k["whole"] else sci(k["ref_err"]), k["own"],
            sci(w["rk4"]), sci(w["euler"]) if "euler" in w else "—", sci(w["expm"]) if "expm" in w else "—"))
    out += [""]
    for k in kern:
        out += ["### %s" % k["name"], ""]
        if k["desc"]:
            out += ["> %s" % k["desc"], ""]
        out += ["![%s: kernel against %s](plots/%s.svg)" % (k["name"], k["ref"], k["name"]), ""]
        if not k["whole"]:
            out += ["Not compared with the Giannantoni engine: the projection carries %s." % k["cov"]]
            if k["reasons"]:
                out += ["What it cannot carry:", ""] + ["- %s" % r for r in k["reasons"]]
            out += [""]
        if len(k["stores"]) > len(k["shown"]):
            out += ["Plotted: %s. Also in the difference: %s." % (
                ", ".join("`%s`" % s for s in k["shown"]),
                ", ".join("`%s`" % s for s in k["stores"][len(k["shown"]):])), ""]
    out += ["## Giannantoni-format models", "",
            "These are in the engine's own format (`examples/giannantoni/`), which the kernel does not read, "
            "so the comparison here is within the engine: the classical network trajectory beside the "
            "incipient output-projection drift. Run with `--steps %d`." % GIA_STEPS, "",
            "| Model | Component | Q(0) | Q(t_end) | max \\|drift\\| |", "|---|---|---:|---:|---:|"]
    for g in gia:
        if "error" in g:
            out.append("| `%s` | did not run: %s | | | |" % (g["name"], " ".join(g["error"])))
            continue
        for c, q0, q1, d in g["rows"]:
            out.append("| `%s` | `%s` | %s | %s | %s |" % (g["name"], c, fmt(q0), fmt(q1), sci(d) if d else "0"))
    out += [""]
    for g in gia:
        if "error" in g:
            continue
        out += ["### %s" % g["name"], "", "%s." % g["title"], "",
                "![%s: Giannantoni engine](plots/gia_%s.svg)" % (g["name"], g["name"]), ""]
    return "\n".join(out) + "\n"


def generate(tmp):
    os.makedirs(os.path.join(tmp, "plots"), exist_ok=True)
    kern = []
    for m in sorted(glob.glob(os.path.join(ROOT, "examples", "*.json"))):
        if os.path.basename(m) == "invalid_model.json":
            continue
        r = kernel_model(m, tmp)
        if r:
            kern.append(r)
    gia = [gia_model(m, tmp) for m in sorted(glob.glob(os.path.join(ROOT, "examples", "giannantoni", "*.json")))]
    with open(os.path.join(tmp, "engine_comparison.md"), "w") as f:
        f.write(page(kern, gia))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="fail if the committed page or plots would change")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        generate(tmp)
        fresh = {os.path.join("plots", p): os.path.join(tmp, "plots", p) for p in os.listdir(os.path.join(tmp, "plots"))}
        fresh["engine_comparison.md"] = os.path.join(tmp, "engine_comparison.md")
        dest = os.path.dirname(OUT_MD)
        if a.check:
            stale = []
            have = {os.path.join("plots", p) for p in os.listdir(OUT_DIR)} if os.path.isdir(OUT_DIR) else set()
            for rel, src in sorted(fresh.items()):
                p = os.path.join(dest, rel)
                if not os.path.exists(p) or open(p).read() != open(src).read():
                    stale.append(rel)
            stale += sorted(have - set(fresh))
            if stale:
                print("check-engine-comparison: FAILED -- stale or missing: %s. Run `make engine-comparison`."
                      % ", ".join(stale))
                return 1
            print("check-engine-comparison: OK -- %d files match a fresh run" % len(fresh))
            return 0
        os.makedirs(OUT_DIR, exist_ok=True)
        for p in os.listdir(OUT_DIR):
            if p.endswith(".svg") and os.path.join("plots", p) not in fresh:
                os.remove(os.path.join(OUT_DIR, p))
        for rel, src in fresh.items():
            with open(src) as s, open(os.path.join(dest, rel), "w") as d:
                d.write(s.read())
        print("wrote %s and %d plots" % (os.path.relpath(OUT_MD, ROOT), len(fresh) - 1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
