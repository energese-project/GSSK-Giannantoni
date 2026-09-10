# GSSK-Giannantoni

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22339312.svg)](https://doi.org/10.5281/zenodo.22339312)
[![CI](https://github.com/energese-project/GSSK-Giannantoni/actions/workflows/deploy.yml/badge.svg)](https://github.com/energese-project/GSSK-Giannantoni/actions/workflows/deploy.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![Language: C99](https://img.shields.io/badge/language-C99-00599C.svg)](include/gssk.h)
[![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20WebAssembly-lightgrey.svg)](#building)
[![Coverage gate](https://img.shields.io/badge/coverage%20gate-%E2%89%A585%25-brightgreen.svg)](#testing)

Two C99 engines over Howard T. Odum's Energy Systems Language: the GSSK kernel, and a second engine implementing Corrado Giannantoni's Incipient Differential Calculus and Maximum Ordinality Principle.

Forked from [energese-project/GSSK](https://github.com/energese-project/GSSK) to pursue the Giannantoni framework independently. The C surface is carried over; the Swift, Python, JavaScript and web trees are not.

A model is a JSON description of storages, sources, sinks and the pathways between them. GSSK integrates the resulting system of ordinary differential equations — Euler, RK4, or adaptive Dormand–Prince — and tracks several carriers independently over the same network: energy, material, money, information. Odum's symbol vocabulary is implemented as node primitives and composite archetypes, with emergy and transformity accounting over the same topology, and deterministic snapshot/replay for reproducible runs. The C API is the surface; WebAssembly is built from the same sources.

`bin/giannantoni_sim` is the second engine. It solves a network by matrix exponential rather than by stepping — `Q(t) = exp(A t) Q(0)`, exact wherever the flow matrix is constant — reports the derivative drift between the incipient and traditional calculi, carries emergy under Odum's non-conservative algebra, and can add components to a graph that is below maximum ordinality. It deliberately shares no headers with the kernel; [ADR 0011](docs/adr/0011-two-engines-declared-lossy-projection.md) records why.

## Quick Demo

```bash
make demo
```

Builds the CLI (if needed) and runs two models — a simple exponential decay and a 4-carrier household ecological-economy — printing the CSV header and first rows of each:

```
=== Decay model (exponential decay, RK4) ===
time,biomass,environment
0.0000,100.000000,0.000000
0.5000,97.530991,2.469009
1.0000,95.122942,4.877058
...

=== Household model (4-carrier ecological-economy) ===
time,salary,bank_account,super_fund,...
0.0000,1.000000,5000.000000,50000.000000,...
... (241 data rows, 24 columns)
```

The decay model follows Q(t) = 100·exp(−0.05·t). The household model has 23 state nodes across money, energy, material, and information carriers.

### Run all regression tests

```bash
make test
```

### The Giannantoni engine

```bash
make demo-giannantoni
```

Runs the same binary three times with no configuration change between them. A seed below maximum ordinality gains a component and reports a generative run; one already closed reports a functional run; feeding the first run's own output back in is a fixed point. Which mode a run was in is decided by diffing the output graph against the seed, not by trusting a flag.

### Benchmark

```bash
make bench
```

---

## Repository Structure

- `include/` — Public API headers (`gssk.h`)
- `src/` — Core C99 implementation
- `bin/` — Compiled executables (CLI tool)
- `lib/` — Compiled libraries (static/shared)
- `tests/` — Regression suite, fuzz target, fuzz corpus
- `examples/` — Reference JSON models
- `bench/` — Benchmark runner and generated models
- `docs/` — Markdown documentation, read here in the repository
- `scripts/` — Release tooling

## Building

### Prerequisites
- GCC or Clang
- Make
- Python 3 (for the schema validator and bench-gen)

### Build

```bash
make            # native library + CLI
make shared     # shared library
make wasm       # WebAssembly (requires emscripten)
```

## Testing

```bash
make test              # regression suite (CSV diff)
make test-giannantoni  # the Giannantoni engine
make bench-giannantoni # incipient closed form vs RK4 and Euler
make test-asan         # AddressSanitizer + UBSan (requires clang)
make coverage-check    # lcov coverage gate ≥ 85%
make test-valgrind     # Valgrind leak check (Linux)
make fuzz-run          # 30 s LibFuzzer run (requires clang)
```

## Documentation

Nothing is published; the documentation is Markdown in this repository and GitHub renders it.

- [docs/concepts.md](docs/concepts.md) — ESL, integration methods, carriers, sensitivity
- [docs/odum_1972_conformance.md](docs/odum_1972_conformance.md) — the two engines scored against Odum's 1972 chapter, module by module, and whether the two frameworks support each other
- [docs/giannantoni_assessment.md](docs/giannantoni_assessment.md) — what the IDC and MOP papers claim, and which claims hold
- [docs/adr/](docs/adr/) — architecture decisions
- [docs/api-reference.md](docs/api-reference.md) — C API
- [docs/cookbook.md](docs/cookbook.md) — parametric sweep, sensitivity, snapshot round-trip
- [docs/CHANGELOG.md](docs/CHANGELOG.md) — release history

## NPM Installation (GitHub)

```bash
npm install energese-project/GSSK-Giannantoni#dist
```

Installs the pre-compiled WASM binaries and TypeScript definitions.

## Releasing

```bash
./scripts/release.sh 3.6.0
git push origin main --tags
```

The script bumps `GSK_VERSION_*` in `include/gssk.h`, updates `docs/CHANGELOG.md`, commits and tags. The GitHub Action builds WASM, updates the `dist` branch, and creates a GitHub Release.

## Citation

If you use GSSK in published work, please cite the archived release. Citation metadata for this repository is in [CITATION.cff](CITATION.cff); GitHub renders it under **Cite this repository**, and Zenodo reads it when minting each deposit.

The badge above resolves to the **concept DOI** — [10.5281/zenodo.22339312](https://doi.org/10.5281/zenodo.22339312) — which always redirects to the most recent release. Cite that when you mean "GSSK" as an ongoing work. To pin a result to the exact code that produced it, cite the **version DOI** instead: every release is minted its own, and they are listed under *Versions* on the [Zenodo record](https://doi.org/10.5281/zenodo.22339312). Reproducibility claims should use the version DOI, because the concept DOI moves with each release.

```bibtex
@software{maud_gssk,
  author    = {Maud, Sholto},
  title     = {{GSSK} --- General Systems Simulation Kernel},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22339312},
  url       = {https://doi.org/10.5281/zenodo.22339312}
}
```

The companion whitepaper — *Three Accountings, One Model: Reconciling Odum's Dynamics, Emergy, and IFRS/AASB Financial Reporting* ([doco/conformance.tex](doco/conformance.tex), and [ADR 0009](docs/adr/0009-accounting-standard-conformance.md)) — states the conformance argument and the boundary conditions on it, including where the kernel is a faithful Odum simulator and where it is not yet a faithful implementation of the emergy algebra. Cite it alongside the software if you rely on that argument.

## Security

See [SECURITY.md](SECURITY.md) for vulnerability reporting and the threat model.

## License

MIT — see [LICENSE](LICENSE). Bundles [cJSON](https://github.com/DaveGamble/cJSON) (MIT).
