// The regression corpus, run through the WASM build.
//
// `make test` runs every model in examples/ through the native CLI and
// compares with tests/expected. This does the same through dist/gssk.wasm
// and the ES module loader, with the same loop as src/main.c and the same
// 1e-6 tolerance as tests/csv_compare.c. Before this, only the forcing
// evaluator was checked under WASM; a kernel change could break the shipped
// artefact for every other model and nothing would fail.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile, readdir } from 'node:fs/promises';

import createGSSK from '../../dist/gssk.mjs';

const root = new URL('../../', import.meta.url);
const TOLERANCE = 1e-6; // tests/csv_compare.c

/** Run a model exactly as src/main.c's `run` does; return the CSV rows as numbers. */
function run(mod, json) {
  const n = mod.lengthBytesUTF8(json) + 1;
  const jp = mod._malloc(n);
  mod.stringToUTF8(json, jp, n);
  const out = mod._malloc(4);
  const status = mod._GSSK_Init(jp, out);
  mod._free(jp);
  const inst = mod.HEAPU32[out >> 2];
  mod._free(out);
  assert.equal(status, 0, 'GSSK_Init failed');

  const size = mod._GSSK_GetStateSize(inst);
  const header = ['time'];
  for (let i = 0; i < size; i++) header.push(mod.UTF8ToString(mod._GSSK_GetNodeID(inst, i)) || 'unknown');

  const rows = [];
  let t = mod._GSSK_GetTStart(inst);
  const tEnd = mod._GSSK_GetTEnd(inst);
  const dt = mod._GSSK_GetDt(inst);
  while (t <= tEnd + dt * 0.01) {
    const state = mod._GSSK_GetState(inst) >> 3;
    rows.push([t, ...mod.HEAPF64.subarray(state, state + size)]);
    if (mod._GSSK_Step(inst, dt) !== 0) break;   // main.c: divergence ends the run
    t += dt;
  }
  mod._GSSK_Free(inst);
  return { header, rows };
}

function parseCsv(text) {
  const [head, ...lines] = text.trim().split('\n');
  return { header: head.split(','), rows: lines.map(l => l.split(',').map(Number)) };
}

const examples = (await readdir(new URL('examples/', root))).filter(f => f.endsWith('.json')).sort();
const mod = await createGSSK();

for (const file of examples) {
  const name = file.slice(0, -'.json'.length);
  let expectedText;
  try {
    expectedText = await readFile(new URL(`tests/expected/${name}.csv`, root), 'utf8');
  } catch {
    continue; // no expected output: the native suite skips it too
  }

  test(`${name} matches tests/expected under WASM`, async () => {
    const json = await readFile(new URL(`examples/${file}`, root), 'utf8');
    const got = run(mod, json);
    const want = parseCsv(expectedText);

    assert.deepEqual(got.header, want.header, 'columns differ');
    assert.equal(got.rows.length, want.rows.length, 'row count differs');
    for (let r = 0; r < want.rows.length; r++) {
      for (let c = 0; c < want.rows[r].length; c++) {
        // The CSV holds %.4f time and %.6f state, so compare at that precision.
        const g = Number(got.rows[r][c].toFixed(c === 0 ? 4 : 6));
        const w = want.rows[r][c];
        assert.ok(Math.abs(g - w) <= TOLERANCE,
          `row ${r + 1}, ${want.header[c]}: WASM ${g} vs expected ${w}`);
      }
    }
  });
}
