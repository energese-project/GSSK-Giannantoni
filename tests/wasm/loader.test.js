// The ES module loader for gssk.wasm: what a consumer can rely on.
//
// The loader replaces Emscripten's generated glue. It keeps the shape
// consumers already use — `_GSSK_*` functions, `_malloc`/`_free`, heap views
// and three string helpers — so the contract is tested here directly rather
// than assumed from a toolchain.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

import createGSSK from '../../dist/gssk.js';

const root = new URL('../../', import.meta.url);

/** The export list, read from the Makefile: one list, so the build and this test cannot disagree. */
async function declaredExports() {
  const mk = await readFile(new URL('Makefile', root), 'utf8');
  const block = mk.slice(mk.indexOf('WASM_EXPORTS'), mk.indexOf('\n\n', mk.indexOf('WASM_EXPORTS')));
  return [...new Set(block.match(/\b(?:GSSK_\w+|malloc|free)\b/g))];
}

test('every export declared in the Makefile is a function on the module', async () => {
  const mod = await createGSSK();
  const names = await declaredExports();
  assert.ok(names.length > 90, `expected the full export list, found ${names.length}`);
  for (const name of names) {
    assert.equal(typeof mod[`_${name}`], 'function', `_${name} missing from the module`);
  }
});

test('a host import the loader does not provide is refused at load, by name', async () => {
  // A minimal module importing env.mystery : () -> ().
  const bytes = new Uint8Array([
    0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00,           // magic, version
    0x01, 0x04, 0x01, 0x60, 0x00, 0x00,                       // type: func () -> ()
    0x02, 0x0f, 0x01, 0x03, 0x65, 0x6e, 0x76,                 // import "env"
    0x07, 0x6d, 0x79, 0x73, 0x74, 0x65, 0x72, 0x79, 0x00, 0x00, // "mystery", func type 0
  ]);
  await assert.rejects(createGSSK({ wasm: bytes }), /env\.mystery/);
});

test('strings cross the boundary intact, including non-ASCII', async () => {
  const mod = await createGSSK();
  const version = mod.UTF8ToString(mod._GSSK_GetVersionString());
  assert.match(version, /^\d+\.\d+\.\d+/);

  const s = 'Δt — Odum’s “emergy”';
  const n = mod.lengthBytesUTF8(s) + 1;
  assert.equal(n, new TextEncoder().encode(s).length + 1);
  const p = mod._malloc(n);
  mod.stringToUTF8(s, p, n);
  assert.equal(mod.UTF8ToString(p), s);
  mod._free(p);
});

test('stringToUTF8 never writes past maxBytes, and always terminates', async () => {
  const mod = await createGSSK();
  const p = mod._malloc(8);
  mod.HEAPU8.fill(0xAA, p, p + 8);
  mod.stringToUTF8('abcdefghij', p, 5);
  assert.equal(mod.UTF8ToString(p), 'abcd');
  assert.equal(mod.HEAPU8[p + 5], 0xAA, 'wrote beyond maxBytes');
  mod._free(p);
});

test('heap views follow the memory when it grows', async () => {
  const mod = await createGSSK();
  const before = mod.HEAPU8.length;
  const big = mod._malloc(64 * 1024 * 1024);
  assert.notEqual(big, 0, 'allocation failed');
  assert.ok(mod.HEAPU8.length > before, 'memory did not grow; the test proves nothing');
  mod.HEAPF64[(big >> 3) + 1000] = 42.5;
  assert.equal(mod.HEAPF64[(big >> 3) + 1000], 42.5);
  mod.HEAPU32[big >> 2] = 0xDEADBEEF;
  assert.equal(mod.HEAPU32[big >> 2], 0xDEADBEEF);
  mod._free(big);
});

test('a model initialises, steps, and reports its nodes through the module', async () => {
  const mod = await createGSSK();
  const json = await readFile(new URL('examples/decay_model.json', root), 'utf8');
  const n = mod.lengthBytesUTF8(json) + 1;
  const jp = mod._malloc(n);
  mod.stringToUTF8(json, jp, n);
  const out = mod._malloc(4);
  assert.equal(mod._GSSK_Init(jp, out), 0);
  const inst = mod.HEAPU32[out >> 2];
  assert.equal(mod._GSSK_GetStateSize(inst), 2);
  assert.equal(mod.UTF8ToString(mod._GSSK_GetNodeID(inst, 0)), 'biomass');
  assert.equal(mod._GSSK_Step(inst, mod._GSSK_GetDt(inst)), 0);
  const q = mod.HEAPF64[mod._GSSK_GetState(inst) >> 3];
  assert.ok(q < 100 && q > 90, `biomass after one step: ${q}`);
  mod._GSSK_Free(inst);
  mod._free(out);
  mod._free(jp);
});
