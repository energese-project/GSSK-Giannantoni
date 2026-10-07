// Vite for `make dev`: serves dev/ and lets it import ../dist/gssk.js and
// ../examples/*.json from the repository root.

import { defineConfig } from 'vite';

const here = new URL('.', import.meta.url).pathname;
const root = new URL('..', import.meta.url).pathname;

export default defineConfig({
  root: here,
  publicDir: false,
  server: { fs: { allow: [root] } },
});
