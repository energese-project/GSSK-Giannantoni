# GSSK-Giannantoni

<!-- ZENODO: disabled until this repository has its own Zenodo record. This DOI is energese-project/GSSK's.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22339312.svg)](https://doi.org/10.5281/zenodo.22339312)
-->
[![CI](https://github.com/energese-project/GSSK-Giannantoni/actions/workflows/deploy.yml/badge.svg)](https://github.com/energese-project/GSSK-Giannantoni/actions/workflows/deploy.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
[![Language: C99](https://img.shields.io/badge/language-C99-00599C.svg)](include/gssk.h)
[![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20WebAssembly-lightgrey.svg)](#building)
[![Coverage gate](https://img.shields.io/badge/coverage%20gate-kernel%20%E2%89%A535%25%20%C2%B7%20Giannantoni%20%E2%89%A590%25-brightgreen.svg)](#testing)

Two C99 engines over Howard T. Odum's Energy Systems Language: the GSSK kernel, and a second engine built to pursue Corrado Giannantoni's Incipient Differential Calculus (IDC) and Maximum Ordinality Principle (MOP).

> **Status.** The second engine implements what [docs/requirements/](docs/requirements/README.md) marks `implemented`, each verified against the sources (ADR 0018):
>
> - the incipient calculus for single-variable equations (`idc.h`);
> - Odum's emergy algebra and Giannantoni's ordinal forms;
> - the MOP's First Fundamental Equation, the printed solution of the Second, and the EQS in the Relational Space algebra (`mop.h`, `relational.h`);
> - Ordinality and Maximum Ordinality as ADR 0021 defines them, with the generative step under the Maximum Em-Power Principle;
> - a harmony detector that decides, per construction, whether harmony is imposed, transported or absent.
>
> Its network **trajectories remain the classical matrix exponential**, not an incipient solution. The harmony *constructor* is labelled `assumed`. Where the sources leave a point open, the reading taken is recorded in PLANLOG.md and pinned by a test. Every output is labelled `implemented`, `classical`, `assumed`, `proxy` or `illustrative` in its run report; [docs/giannantoni_assessment.md](docs/giannantoni_assessment.md#status-of-record-planmd-2) is the status of record, and PLAN.md §6 lists what the sources do not let it compute.

Forked from [energese-project/GSSK](https://github.com/energese-project/GSSK) to pursue the Giannantoni framework independently. The C surface is carried over; the Swift, Python, JavaScript and web trees are not.

A model is a JSON description of storages, sources, sinks and the pathways between them. GSSK integrates the resulting system of ordinary differential equations — Euler, RK4, or adaptive Dormand–Prince — and tracks several carriers independently over the same network: energy, material, money, information. Odum's symbol vocabulary is implemented as node primitives and composite archetypes, with emergy and transformity accounting over the same topology, and deterministic snapshot/replay for reproducible runs. The C API is the surface; WebAssembly is built from the same sources.

`bin/giannantoni_sim` is the second engine. It solves a network by matrix exponential rather than by stepping — `Q(t) = exp(A t) Q(0)`, exact wherever the flow matrix is constant, and classical rather than incipient — carries emergy under Odum's non-conservative algebra (implemented), reports solution drift where the sources define it (exactly zero for a constant flow matrix, refused otherwise) and the output-projection drift of each trajectory (implemented), and, below Maximum Ordinality (strong connectivity of the component graph), adds pathways one at a time, choosing each by maximum total empower (ADR 0021). It deliberately shares no headers with the kernel; [ADR 0011](docs/adr/0011-two-engines-declared-lossy-projection.md) records why.

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
... (241 data rows, 25 columns)
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

Runs the same binary three times with no configuration change between them. A seed below Maximum Ordinality gains pathways and reports a generative run; one already closed reports a functional run; feeding the first run's own output back in is a fixed point. Which mode a run was in is decided by diffing the output graph against the seed, not by trusting a flag.

Five seeds live in `examples/giannantoni/`, each showing something the others do not:

| Seed | Shows |
|---|---|
| `input.json` | Odum's work gate. Below maximum, one pathway added. |
| `closed_loop.json` | Already at Maximum Ordinality: a functional run. |
| `trophic_chain.json` | A network harmony verdict actually computed: the pathways carry emergy, and its `mop` block sets the reference couple. |
| `coproduction.json` | ½ and 2 couples in the ordinality record. Three pathways are added. |
| `harmonic_couples.json` | A `mop` block whose First Equation row is harmonic: `--mop-out` writes `R_H = 0`. |

Each seed's verdicts are in [docs/results/harmony_verdicts.md](docs/results/harmony_verdicts.md), and `tests/mop_cli.sh` holds the CLI to them.

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
make wasm       # WebAssembly: dist/gssk.wasm + dist/gssk.js (fetches the pinned WASI SDK)
```

### In the browser, locally

```bash
make dev        # http://localhost:5173/ — any example model, run through gssk.js + gssk.wasm
```

Needs only the Apple `container` CLI. `gssk.wasm` is built natively by `make wasm` (the pinned WASI SDK) and Vite runs in a Node container. To run the page on the bytes a release would ship instead, `make dev-guix` takes Node and `gssk.wasm` from the pinned Guix (`guix/`); its first run fetches Guix packages (about 15 minutes), and `make guix-down` stops that container and keeps its store.

## Testing

```bash
make test              # regression suite (CSV diff)
make test-giannantoni  # the Giannantoni engine
make bench-giannantoni # incipient closed form vs RK4 and Euler
make test-asan         # AddressSanitizer + UBSan (requires clang)
make coverage-check    # line-coverage gates: kernel ≥ 35% (lcov), Giannantoni units ≥ 90% (gcov)
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

If you use GSSK in published work, please cite it. Citation metadata for this repository is in [CITATION.cff](CITATION.cff); GitHub renders it under **Cite this repository**.

<!-- ZENODO: disabled until this repository has its own Zenodo record; the DOI
below is energese-project/GSSK's. Restore with the new record's DOI.

Zenodo reads CITATION.cff when minting each deposit.

The badge above resolves to the **concept DOI** — [10.5281/zenodo.22339312](https://doi.org/10.5281/zenodo.22339312) — which always redirects to the most recent release. Cite that when you mean "GSSK" as an ongoing work. To pin a result to the exact code that produced it, cite the **version DOI** instead: every release is minted its own, and they are listed under *Versions* on the [Zenodo record](https://doi.org/10.5281/zenodo.22339312). Reproducibility claims should use the version DOI, because the concept DOI moves with each release.

  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22339312},
  url       = {https://doi.org/10.5281/zenodo.22339312}
-->

```bibtex
@software{maud_gssk,
  author    = {Maud, Sholto},
  title     = {{GSSK} --- General Systems Simulation Kernel},
  year      = {2026},
  url       = {https://github.com/energese-project/GSSK-Giannantoni}
}
```

The companion whitepaper — *Three Accountings, One Model: Reconciling Odum's Dynamics, Emergy, and IFRS/AASB Financial Reporting* ([doco/conformance.tex](doco/conformance.tex), and [ADR 0009](docs/adr/0009-accounting-standard-conformance.md)) — states the conformance argument and the boundary conditions on it, including where the kernel is a faithful Odum simulator and where it is not yet a faithful implementation of the emergy algebra. Cite it alongside the software if you rely on that argument.

## Security

See [SECURITY.md](SECURITY.md) for vulnerability reporting and the threat model.

## License

MIT — see [LICENSE](LICENSE). Bundles [cJSON](https://github.com/DaveGamble/cJSON) (MIT).
