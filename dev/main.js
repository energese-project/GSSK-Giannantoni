// The dev page: run any example model through dist/gssk.js in the browser.
//
// This is the loader's browser path — fetch, streaming compilation, the
// application/wasm type — which the Node test suites cannot exercise.

import '../docs/.vitepress/theme/energese-tokens.css';
import createGSSK from '../dist/gssk.js';

const models = import.meta.glob('../examples/*.json', { query: '?raw', import: 'default' });
const $ = (id) => document.getElementById(id);

// The series palette is validated in order for colour-blind safety (AGENTS.md):
// past eight series, hue repeats and line style carries identity instead.
const DASHES = ['', '7 4', '2 3', '9 3 2 3'];

/** Run a model as src/main.c and tests/wasm/corpus.test.js do. */
function run(gssk, json) {
  const n = gssk.lengthBytesUTF8(json) + 1;
  const jp = gssk._malloc(n);
  gssk.stringToUTF8(json, jp, n);
  const out = gssk._malloc(4);
  const status = gssk._GSSK_Init(jp, out);
  gssk._free(jp);
  const inst = gssk.HEAPU32[out >> 2];
  gssk._free(out);
  if (status !== 0) {
    const why = inst ? gssk.UTF8ToString(gssk._GSSK_GetErrorDescription(inst)) : `status ${status}`;
    if (inst) gssk._GSSK_Free(inst);
    throw new Error(`GSSK_Init failed: ${why}`);
  }
  const size = gssk._GSSK_GetStateSize(inst);
  const names = Array.from({ length: size }, (_, i) => gssk.UTF8ToString(gssk._GSSK_GetNodeID(inst, i)));
  const rows = [];
  let t = gssk._GSSK_GetTStart(inst);
  const tEnd = gssk._GSSK_GetTEnd(inst);
  const dt = gssk._GSSK_GetDt(inst);
  let note = '';
  while (t <= tEnd + dt * 0.01) {
    const p = gssk._GSSK_GetState(inst) >> 3;
    rows.push([t, ...gssk.HEAPF64.subarray(p, p + size)]);
    if (gssk._GSSK_Step(inst, dt) !== 0) { note = `stopped at t=${t.toFixed(4)}: numerical divergence`; break; }
    t += dt;
  }
  gssk._GSSK_Free(inst);
  return { names, rows, note };
}

function chart({ names, rows }) {
  const W = 900, H = 360, L = 56, R = 12, T = 12, B = 32;
  const xs = rows.map(r => r[0]);
  const ys = rows.flatMap(r => r.slice(1)).filter(Number.isFinite);
  const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = [Math.min(0, ...ys), Math.max(...ys)];
  if (y0 === y1) y1 = y0 + 1;
  const X = x => L + (x - x0) / (x1 - x0 || 1) * (W - L - R);
  const Y = y => H - B - (y - y0) / (y1 - y0) * (H - T - B);
  const css = getComputedStyle(document.documentElement);
  const colour = i => css.getPropertyValue(`--e-series-${(i % 8) + 1}`).trim();
  const dash = i => DASHES[Math.floor(i / 8) % DASHES.length];

  const ticks = [0, 0.25, 0.5, 0.75, 1];
  let svg = '';
  for (const f of ticks) {
    const y = y0 + f * (y1 - y0), x = x0 + f * (x1 - x0);
    svg += `<line x1="${L}" x2="${W - R}" y1="${Y(y)}" y2="${Y(y)}" stroke="var(--e-grid)"/>`;
    svg += `<text x="${L - 6}" y="${Y(y) + 4}" text-anchor="end" font-size="11" fill="var(--e-tick)">${+y.toPrecision(4)}</text>`;
    svg += `<text x="${X(x)}" y="${H - B + 16}" text-anchor="middle" font-size="11" fill="var(--e-tick)">${+x.toPrecision(4)}</text>`;
  }
  svg += `<line x1="${L}" x2="${W - R}" y1="${H - B}" y2="${H - B}" stroke="var(--e-axis)"/>`;
  names.forEach((_, i) => {
    const d = rows.map((r, k) => `${k ? 'L' : 'M'}${X(r[0]).toFixed(1)},${Y(r[i + 1]).toFixed(1)}`).join('');
    svg += `<path d="${d}" fill="none" stroke="${colour(i)}" stroke-width="1.75" stroke-dasharray="${dash(i)}"/>`;
  });
  $('chart').innerHTML = svg;

  // A legend is always present: identity must not rest on colour alone.
  $('legend').innerHTML = names.map((name, i) =>
    `<span><svg width="24" height="8" aria-hidden="true"><line x1="0" x2="24" y1="4" y2="4" stroke="${colour(i)}" stroke-width="2" stroke-dasharray="${dash(i)}"/></svg>${name}</span>`
  ).join('');
}

function table({ names, rows }) {
  const head = `<tr><th>time</th>${names.map(n => `<th>${n}</th>`).join('')}</tr>`;
  const body = rows.map(r => `<tr>${r.map((v, c) => `<td>${v.toFixed(c ? 6 : 4)}</td>`).join('')}</tr>`).join('');
  $('table').innerHTML = head + body;
}

async function main() {
  const status = $('status');
  for (const path of Object.keys(models).sort()) {
    const opt = document.createElement('option');
    opt.value = path;
    opt.textContent = path.split('/').pop().replace(/\.json$/, '');
    $('model').append(opt);
  }
  $('model').value = Object.keys(models).find(p => p.endsWith('/decay_model.json')) ?? $('model').value;

  const gssk = await createGSSK({ printErr: line => console.warn('[gssk]', line) });
  $('kernel').textContent = `Kernel ${gssk.UTF8ToString(gssk._GSSK_GetVersionString())}.`;

  const go = async () => {
    status.className = '';
    status.textContent = 'running…';
    try {
      const json = await models[$('model').value]();
      const t0 = performance.now();
      const result = run(gssk, json);
      const ms = (performance.now() - t0).toFixed(1);
      chart(result);
      table(result);
      status.textContent = `${result.rows.length} steps in ${ms} ms${result.note ? ` — ${result.note}` : ''}`;
    } catch (err) {
      status.className = 'error';
      status.textContent = err.message;
    }
  };
  $('run').addEventListener('click', go);
  $('model').addEventListener('change', go);
  go();
}

main().catch(err => {
  $('status').className = 'error';
  $('status').textContent = `Could not load gssk.wasm: ${err.message}`;
});
