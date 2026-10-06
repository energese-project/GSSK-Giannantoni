# Spike: a reproducible GSSK WASM build with Guix

**Question.** Can Guix replace the container images as the pinned, archivable build
environment for GSSK's WebAssembly kernel and the TypeScript tools around it, and does
it solve the "Apple `container` vs Docker vs Podman" problem?

**Short answer.** Yes for the C/WASM kernel, with Emscripten dropped. Mostly for
TypeScript, with one version gap. Guix does not remove the need for a container runtime
on macOS, but it removes the runtime as a source of difference: the image is built by
Guix, not by the runtime.

## What is here

| File | Purpose |
|---|---|
| `gssk.scm` | Guix packages: wasm32 compiler-rt builtins, wasi-libc, and `gssk-wasm` (library + CLI). The `check` phase runs every regression model under wasm3 |
| `channels.scm` | pins Guix to commit `230aa37` (release 1.5.0) for `guix time-machine` |
| `wasm.sha256` | the WASM hashes produced on aarch64; CI on x86_64 must match them |
| `Containerfile.guix-host`, `entrypoint.sh` | a plain Debian box that only hosts Guix, for macOS |
| `../../.github/workflows/guix-spike.yml` | the x86_64 rebuild, `--check`, and Docker run |

```sh
guix build -f spike/guix/gssk.scm                                   # build
guix time-machine -C spike/guix/channels.scm -- build -f spike/guix/gssk.scm   # pinned
```

## Results

| Test | Result |
|---|---|
| Guix toolchain available prebuilt for aarch64 (clang 21.1.5, lld, llvm, node 22.14, esbuild, wabt, wasm3, skopeo) | yes, 100% substitutes |
| Emscripten in Guix | **no** |
| WASI libc in Guix | **no**, defined here from `wasi-sdk-34` |
| TypeScript compiler in Guix | **no** (esbuild is, and can strip types) |
| GSSK compiles to `wasm32-wasip1` with `-Wall -Wextra -Werror` | yes, no warnings |
| 23 regression models, run as WASM under wasm3, vs `tests/expected` | **23/23 byte-identical** |
| Clean rebuilds of all three packages (delete, rebuild, compare nar hash) | **6 consecutive identical** (see caveat 4) |
| `guix pack` image rebuilt from scratch | identical tarball |
| Same image converted to OCI twice with skopeo | identical manifest and layer digests; outer tar differs by timestamps only |
| Image under Apple `container`: regression suite | **23/23 byte-identical** |
| aarch64 (Apple container) vs x86_64 (GitHub Actions) WASM bytes | **identical** |
| `guix build --check` with build isolation, all three packages (x86_64) | **pass**: each rebuilt bit-identical |
| Guix-packed image under Docker (x86_64): regression suite | **23/23 byte-identical**; Docker loads the Guix pack directly |
| Starting from Ubuntu's Guix 1.4.0, `time-machine` to the pinned 1.5.0 | works; whole CI job 34 minutes |
| Dashboard (TypeScript) on Guix's Node 22: `npm ci`, `tsc`, `vite build` | pass |
| Dashboard tests on Node 22 with `--experimental-strip-types` | 215/218; the 3 failures read `.github/ISSUE_TEMPLATE/add-program.yml`, which the GSSK-Dashboard copy lacks. Not a Node problem |

The kernel needs only four host imports (`clock_time_get`, `fd_close`, `fd_seek`,
`fd_write`), so a browser host is a few lines of JavaScript, with no WASI runtime.

## Findings that shape the decision

1. **Dropping Emscripten is cheap for the kernel and costs a loader.** The build is plain
   clang. But `createGSSK`, `ccall`, `cwrap` and the heap views are Emscripten's JS glue,
   used by `src/gssk.d.ts`, `tests/wasm/forcing_parity.cjs` and the API docs. A small
   hand-written loader replaces them. It is the one real piece of migration work.
2. **Guix lags upstream for JavaScript.** Node is 22.14; the dashboard declares
   `node >= 26`. TypeScript is not packaged, and npm dependencies come from the npm
   registry pinned by `package-lock.json` integrity hashes, not built by Guix. So the
   TS tools are *hash-pinned*, not *bootstrappable*. Archiving the built JS covers re-running;
   rebuilding from source depends on the npm registry.
3. **On macOS, Guix still runs inside a container, and Apple's blocks build isolation.**
   The daemon needs `--disable-chroot` there unless the container is granted extra
   capabilities. Without isolation, `guix build --check` cannot run, `/tmp` must be
   world-writable, and local builds are less hermetic. Substitutes and their hashes are
   unaffected. A Linux machine or CI runs the daemon normally.
4. **One unexplained hash.** The first hash of the builtins package, taken immediately
   after a failed `--check` attempt under `--disable-chroot`, differed from the next six,
   which all agree. The GSSK WASM linked from it was identical either way. The isolated
   `--check` in CI passed for all three packages, so the outlier is most likely an artefact
   of the failed non-isolated `--check`, but that is inferred, not shown.
5. **Runtimes differ in archive format, not image format.** `guix pack -f docker` writes
   a Docker archive; Apple's `container image load` accepts only OCI layout. `skopeo`
   (in Guix) converts, keeping the digests. Docker loads the Guix pack directly.
6. **Guix's clang leaks host headers into cross builds.** Its driver adds glibc's include
   path even with `--sysroot`. Every wasm32 compile needs `-nostdlibinc`, and the build
   environment's `C_INCLUDE_PATH` must be unset. `gssk.scm` does both.
7. **Standalone compiler-rt builds are deprecated upstream.** LLVM warns they will become
   an error. A future LLVM bump should take compiler-rt from Guix's monorepo source.
8. **GNU's infrastructure is a single point of failure for *installing* Guix.** One
   `cmake-minimal` substitute download stalled for an hour, and on 2026-10-06
   `ftp.gnu.org` and its mirror redirector were unreachable for the whole session, so
   CI installs Guix from Ubuntu's archive instead. Guix's own git and substitute servers
   stayed up. `--fallback` belongs in every scripted build, and the archive should hold
   the Guix installer and the pack image, not only point at where to download them.

## Recommendation

- **Kernel:** adopt the Guix build, drop Emscripten, write the JS loader. `wasm.sha256`
  becomes a release artefact, and CI fails when a rebuild disagrees with it.
- **Archive per release:** `channels.scm` + `gssk.scm` (rebuild), the WASM files and
  their hashes (re-run), the `guix pack` image (run without Guix), and the regression
  corpus (verify).
- **TypeScript tools:** build with Guix's Node, pin npm by lockfile, archive the built JS.
  Write the tools against Node 22 or move the `engines` floor down, since Guix's Node is
  what the archive will reproduce.
- **Container runtimes:** stop building images with each runtime's builder. Guix builds
  the image; Docker loads it directly, Apple `container` after `skopeo` converts it.
