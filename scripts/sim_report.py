#!/usr/bin/env python3
"""The simulation regression tests, as a GitHub-flavoured markdown report.

`make test` says PASSED or FAILED per model and stops at the first failure.
That is the right gate, but it shows nothing of what the simulations did. This
report runs every model in examples/ through bin/gssk, compares each run with
its golden trajectory using bin/csv_compare (the gate's own comparator, so the
report and the gate cannot disagree), and writes, for every model:

  - the verdict, and the largest deviation from the golden file;
  - each node's start, minimum, maximum and end, with a sparkline of its
    trajectory, so a reviewer can see that a decay decays and a cycle cycles;
  - the model's own description of what it is meant to show.

It does not stop at a failure: every model is reported. It exits 1 if any
model failed, so that a CI step running it after the gate still shows red.

    python3 -I scripts/sim_report.py                 # markdown to stdout
    python3 -I scripts/sim_report.py --expected DIR  # another golden dir (self-test)

In CI: `make -s sim-report >> "$GITHUB_STEP_SUMMARY"`.

Verifies: FR-OUT-004 (T-OUT-04)
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
GSSK = os.path.join(ROOT, "bin", "gssk")
COMPARE = os.path.join(ROOT, "bin", "csv_compare")
TOLERANCE = 1e-6          # tests/csv_compare.c TOLERANCE, quoted in the report
BARS = "▁▂▃▄▅▆▇█"
SPARK_WIDTH = 24


def read_csv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return [], []
    header = rows[0]
    data = [[float(x) for x in r] for r in rows[1:] if r]
    return header, data


def sparkline(values):
    if not values:
        return ""
    n = len(values)
    idx = [round(i * (n - 1) / (SPARK_WIDTH - 1)) for i in range(SPARK_WIDTH)] if n > SPARK_WIDTH else range(n)
    pts = [values[i] for i in idx]
    lo, hi = min(pts), max(pts)
    if hi - lo <= 1e-12 * max(1.0, abs(hi)):
        return "flat"
    return "".join(BARS[min(len(BARS) - 1, int((v - lo) / (hi - lo) * len(BARS)))] for v in pts)


def num(x):
    if x == 0:
        return "0"
    if abs(x) >= 1e5 or abs(x) < 1e-3:
        return f"{x:.3e}"
    return f"{x:.4g}"


def allowlist(path):
    out = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                name, _, why = line.partition("#")
                out[name.strip()] = why.strip()
    except OSError:
        pass
    return out


def max_deviation(expected, got):
    """Largest |expected - got| over every shared cell, or None if the shapes differ."""
    he, de = read_csv(expected)
    hg, dg = read_csv(got)
    if he != hg or len(de) != len(dg):
        return None
    worst = 0.0
    for re_, rg in zip(de, dg):
        if len(re_) != len(rg):
            return None
        for a, b in zip(re_, rg):
            worst = max(worst, abs(a - b))
    return worst


def run_model(model, out_csv):
    r = subprocess.run([GSSK, model, out_csv], capture_output=True, text=True)
    return r.returncode, (r.stderr or r.stdout).strip()


def report(expected_dir, skip_path):
    skips = allowlist(skip_path)
    rows, details, failed = [], [], 0
    counts = {"passed": 0, "failed": 0, "skipped": 0}
    with tempfile.TemporaryDirectory() as tmp:
        for model in sorted(glob.glob(os.path.join(ROOT, "examples", "*.json"))):
            name = os.path.splitext(os.path.basename(model))[0]
            with open(model) as f:
                doc = json.load(f)
            meta = doc.get("metadata", {})
            cfg = doc.get("config", {})
            types = {n.get("id"): n.get("type", "") for n in doc.get("nodes", [])}
            out_csv = os.path.join(tmp, name + ".csv")
            rc, msg = run_model(model, out_csv)
            golden = os.path.join(expected_dir, name + ".csv")

            if name in skips and not os.path.exists(golden):
                counts["skipped"] += 1
                rows.append(f"| `{name}` | — | — | — | ⏭️ skipped: {skips[name]} | |")
                continue
            if rc != 0 or not os.path.exists(out_csv):
                counts["failed"] += 1
                rows.append(f"| `{name}` | — | — | — | ❌ did not run: {msg.splitlines()[0] if msg else rc} | |")
                continue
            if not os.path.exists(golden):
                counts["failed"] += 1
                verdict, dev = "❌ no golden file (`make test-update`, ADR 0018)", "—"
            else:
                ok = subprocess.run([COMPARE, golden, out_csv], capture_output=True).returncode == 0
                d = max_deviation(golden, out_csv)
                dev = "shape differs" if d is None else ("0" if d == 0 else f"{d:.1e}")
                if ok:
                    counts["passed"] += 1
                    verdict = "✅ matches golden"
                else:
                    counts["failed"] += 1
                    verdict = "❌ differs from golden"

            header, data = read_csv(out_csv)
            cols = list(zip(*data)) if data else []
            steps = max(0, len(data) - 1)
            method = cfg.get("method", "rk4")
            # The headline: storages, since a held source only shows its forcing.
            stores = [i for i, h in enumerate(header[1:], 1) if types.get(h) == "storage"] or list(range(1, len(header)))
            head = [f"{header[i]} {sparkline(cols[i])}" for i in stores[:2]] if cols else []
            more = f" +{len(stores) - 2}" if len(stores) > 2 else ""
            rows.append(f"| `{name}` | {method} | {steps} | {dev} | {verdict} | {'<br>'.join(head)}{more} |")

            lines = [f"<details><summary><code>{name}</code> — {verdict}</summary>", ""]
            desc = meta.get("description")
            if desc:
                lines += [f"> {desc}", ""]
            t0, t1 = (cols[0][0], cols[0][-1]) if cols else (0, 0)
            lines += [f"`{method}`, dt {cfg.get('dt', '—')}, t {num(t0)} → {num(t1)}, {steps} steps. "
                      f"Largest deviation from `{os.path.relpath(golden, ROOT)}`: {dev} "
                      f"(the gate's tolerance is {TOLERANCE:g} per value).", "",
                      "| Node | Type | Start | Min | Max | End | Trajectory |",
                      "|---|---|---:|---:|---:|---:|---|"]
            for i in range(1, len(header)):
                c = cols[i]
                lines.append(f"| `{header[i]}` | {types.get(header[i], '—')} | {num(c[0])} | {num(min(c))} "
                             f"| {num(max(c))} | {num(c[-1])} | {sparkline(c)} |")
            lines += ["", "</details>", ""]
            details.append("\n".join(lines))

    total = sum(counts.values())
    icon = "✅" if counts["failed"] == 0 else "❌"
    out = [
        "## Simulation tests",
        "",
        f"{icon} **{counts['passed']}** of {total} models match their golden trajectories, "
        f"**{counts['failed']}** failed, **{counts['skipped']}** skipped (allowlisted in "
        "`tests/skip_allowlist.txt`).",
        "",
        "Each model under `examples/` is run by `bin/gssk` and compared with `tests/expected/` by "
        f"`bin/csv_compare`, the same comparator `make test` gates on (|Δ| ≤ {TOLERANCE:g} per value). "
        "A golden file pins a trajectory; whether that trajectory is *right* is checked by the "
        "hand-computed suites (`test-forcing`, `test-giannantoni`, …) and by `docs/results/level1_survey.md` "
        "against the exact solution. The sparklines show storages over the run, sampled to "
        f"{SPARK_WIDTH} points and scaled to each node's own range.",
        "",
        "| Model | Method | Steps | Max \\|Δ\\| vs golden | Result | Storages |",
        "|---|---|---:|---:|---|---|",
        *rows,
        "",
        "### Per model",
        "",
        *details,
    ]
    return "\n".join(out) + "\n", counts["failed"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--expected", default=os.path.join(ROOT, "tests", "expected"))
    ap.add_argument("--skip-allowlist", default=os.path.join(ROOT, "tests", "skip_allowlist.txt"))
    a = ap.parse_args()
    text, failed = report(a.expected, a.skip_allowlist)
    sys.stdout.write(text)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
