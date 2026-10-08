# Reproducing a release's `gssk.wasm`

Every tagged release's `gssk.wasm` is built by [Guix](https://guix.gnu.org) from
the two files here, and the toolchain that built it is released beside it. So a
release can be **re-run** (load `gssk.wasm` in any WebAssembly runtime) and
**rebuilt** to the same bytes, on any machine, for as long as either Guix or the
archived toolchain can be run.

| File | What it pins |
|---|---|
| `channels.scm` | the exact Guix commit, so every package is the one that built the release |
| `gssk.scm` | wasi-libc (`wasi-sdk-34`) and wasm32 compiler-rt from pinned upstream source, then GSSK through the Makefile's own `make wasm` and `make test-wasm` |

Everyday builds (`make wasm`) use the WASI SDK and produce different bytes; only
the Guix build is the release artefact.

## Rebuilding a release

Check out the release's tag, then use whichever of these still works.

**1. With Guix** (the recipe, rebuilt from source):

```sh
guix time-machine -C guix/channels.scm -- build -f guix/gssk.scm
sha256sum /gnu/store/…-gssk-wasm-release/share/gssk/gssk.wasm
```

**2. With the archived toolchain** (no Guix; any x86_64 Linux):

```sh
mkdir tc && tar xf gssk-toolchain-x86_64-linux.tar.xz -C tc
tc/bin/make wasm-toolchain TC=$PWD/tc
sha256sum dist/gssk.wasm
```

**3. Without either**: `gssk.wasm` is plain C99 compiled for `wasm32-wasip1`.
Any clang with a wasm32 target and a WASI libc builds it from `src/` with the
flags in the Makefile's `wasm` rule. The bytes may differ; the regression
corpus (`make test-wasm`) says whether the behaviour does.

In 1 and 2 the hash must equal the release's `gssk-guix.sha256`.

## What CI checks

`.github/workflows/guix.yml`, on every version tag (where `deploy.yml` then
releases its output), on PRs that touch this directory, monthly, and on demand —
run it by hand for a PR that changes the Makefile's `wasm` rule:

1. The package's check phase runs `make test-wasm`: the loader, every
   regression model, and forcing parity — the suites the everyday CI runs.
2. `guix build --check` rebuilds each package and fails on any differing bit.
3. x86_64 and aarch64 runners must produce byte-identical `gssk.wasm`.

On tags it also packs the toolchain (`guix pack -RR -C xz`, about 450 MB): clang
and lld 21, wasi-libc, compiler-rt, make, Node and coreutils, relocatable so it
runs without Guix.

## How this was established

A spike (PR #23) tested the approach before it was adopted:

- GSSK built with Guix's clang and the two packages above, with no Emscripten,
  passed every regression model byte-identically.
- Six clean rebuilds were bit-identical; `--check` with build isolation passed
  in CI; x86_64 and aarch64 produced the same `gssk.wasm`.
- The packed toolchain, unpacked in a plain Debian container with no Guix,
  rebuilt the identical `gssk.wasm` and passed all 37 WASM tests.

## Things that will matter later

- **On macOS, Guix runs in a container** (`Containerfile.guix-host`), and
  Apple's runtime blocks build isolation, so the daemon there uses
  `--disable-chroot`. Fine for trying the recipe; `--check` needs Linux or CI.
- **Guix's clang adds glibc's headers** to every compile, even with `--sysroot`.
  The wasm32 builds pass `-nostdlibinc` and keep the host's `C_INCLUDE_PATH` out.
- **GCC 14 rejects `src/gssk.c` under `-Werror`** (`-Wformat-truncation`), so the
  recipe's native build uses clang. CI's older GCC does not report it yet.
- **A standalone compiler-rt build is deprecated upstream**; LLVM says it will
  become an error. A future LLVM bump should take compiler-rt from Guix's
  monorepo source instead of the release tarball.
- **GNU's file servers can be unreachable** while Guix's own servers are up.
  CI installs Guix from Ubuntu's archive and moves to the pin with
  `time-machine`, and every scripted build passes `--fallback`.
- **Changing the pin** (`guix describe -f channels > guix/channels.scm` on a newer
  Guix) changes the compiler, so it changes the bytes: a new pin belongs in its
  own PR, and the release after it records new hashes.
