/* WASM parity for forcing — requirement 3 of the upstream proposal.
 *
 * There is ONE distributed artifact, and it must not silently differ from
 * native. sin() and exp() are the two that could: a WASM build could
 * substitute a different libm, and the divergence would be small enough to
 * look like round-off while being systematic.
 *
 * MEASURED, rather than assumed. Comparing the same eight sample points three
 * ways — macOS Apple clang, Linux GCC, and WASM — the WASM build agreed with
 * Linux GCC EXACTLY, bit for bit, on all eight for every waveform. macOS is
 * the outlier: its libm sin() differs by 1 ULP at two of the eight points
 * (t=2.5 and t=19.9 on the sine case, ~8.9e-16 and ~2.2e-16).
 *
 * So the tolerance below is NOT slack for WASM. It is slack for the platform
 * this check happens to run on. sin/exp are not required by IEEE-754 to be
 * correctly rounded, so any two libms may legitimately differ by <1 ULP, and
 * demanding bit-equality would mean GSSK shipping its own transcendentals.
 * That is tracked as deterministic-transcendentals-cross-platform.
 *
 * 4 ULP is comfortably above the 1 ULP observed and far below anything that
 * could hide a substituted implementation, a wrong formula, or a dropped
 * parameter — which are the failures this check is actually for.
 *
 * The native side is tests/results/forcing_native.json, written by
 * bin/dump_forcing_native (`make test-wasm` produces it first).
 *
 * Ported from forcing_parity.cjs when the Emscripten build was replaced by
 * clang + wasi-libc and an ES module loader. The libm is now wasi-libc's
 * (musl), so the measurement above is re-made by this test on every run.
 */

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

import createGSSK from '../../dist/gssk.mjs';

const NATIVE = JSON.parse(await readFile(
  new URL('../results/forcing_native.json', import.meta.url), 'utf8'));

// Distance in representable doubles, so the bound is scale-free — an absolute
// epsilon would be far too loose near 1e6 and far too tight near 1e-6.
const view = new DataView(new ArrayBuffer(8));
function ordinal(x) {
  view.setFloat64(0, x);
  const hi = view.getUint32(0), lo = view.getUint32(4);
  const n = (BigInt(hi) << 32n) | BigInt(lo);
  // Map the sign-magnitude layout onto a monotone integer line.
  return (hi & 0x80000000) ? -(n & 0x7fffffffffffffffn) : n;
}
function ulpsApart(a, b) {
  if (a === b) return 0;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return Infinity;
  const d = ordinal(a) - ordinal(b);
  return Number(d < 0n ? -d : d);
}

const ALLOWANCE = 4;
const mod = await createGSSK();

function init(json) {
  const n = mod.lengthBytesUTF8(json) + 1;
  const jp = mod._malloc(n);
  mod.stringToUTF8(json, jp, n);
  const out = mod._malloc(4);
  const status = mod._GSSK_Init(jp, out);
  const inst = mod.HEAPU32[out >> 2];
  mod._free(out);
  mod._free(jp);
  assert.equal(status, 0, 'WASM init failed');
  return inst;
}

function near(got, want, where) {
  assert.ok(ulpsApart(got, want) <= ALLOWANCE,
    `${where}: WASM ${got} vs native ${want} — ${ulpsApart(got, want)} ULP apart, ` +
    `above the ${ALLOWANCE} ULP libm allowance. That is too far to be rounding: ` +
    `suspect a substituted implementation, a wrong formula, or a parameter that ` +
    `did not reach the evaluator.`);
}

for (const cas of NATIVE.cases) {
  test(`node forcing: ${cas.name}`, () => {
    const inst = init(cas.model);
    assert.equal(mod._GSSK_GetNodeForcingKind(inst, 0), cas.kind, 'forcing kind');
    for (const s of cas.samples) near(mod._GSSK_EvaluateNodeForcing(inst, 0, s.t), s.v, `t=${s.t}`);
    mod._GSSK_Free(inst);
  });
}

// The edge evaluator too, so both attachment points are covered.
test('edge forcing', () => {
  const inst = init(NATIVE.edge_case.model);
  assert.equal(mod._GSSK_GetEdgeForcingKind(inst, 0), NATIVE.edge_case.kind, 'forcing kind');
  for (const s of NATIVE.edge_case.samples) near(mod._GSSK_EvaluateEdgeForcing(inst, 0, s.t), s.v, `t=${s.t}`);
  mod._GSSK_Free(inst);
});
