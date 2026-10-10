# GSSK Build System
CC      ?= gcc
CFLAGS  = -Wall -Wextra -Werror -std=c99 -Iinclude -fPIC
LDFLAGS = -lm

# Architecture flags.
#
# -march=native was previously unconditional. It bakes the build machine's CPU
# into the binary, which is wrong for anything distributed, and it is actively
# fragile: clang rejects some auto-detected feature combinations outright, so a
# CI job landing on an AVX10.1-capable runner fails with
#   error: invalid feature combination: +avx10.1-256 ... [-Winvalid-feature-combination]
# under -Werror, while the identical source builds fine on an older runner.
# That made green-ness depend on which machine picked the job up.
#
# Default is now portable. Opt in with NATIVE=1 for local benchmarking, where
# tuning to the host is the point and reproducibility across machines is not.
ifeq ($(NATIVE), 1)
	ARCH_FLAGS = -march=native
else
	ARCH_FLAGS =
endif

# Optimization levels (Use 'make DEBUG=1' for debugging)
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0 -DDEBUG
else
	CFLAGS += -O3 $(ARCH_FLAGS)
endif

# Directories
SRC_DIR = src
INC_DIR = include
BIN_DIR = bin
LIB_DIR = lib
DIST_DIR = dist
TEST_DIR = tests

# Files
SOURCES = $(SRC_DIR)/gssk.c $(SRC_DIR)/advanced.c $(SRC_DIR)/cJSON.c
OBJECTS = $(LIB_DIR)/gssk.o $(LIB_DIR)/advanced.o $(LIB_DIR)/cJSON.o
TARGET_LIB = $(LIB_DIR)/libgssk.a
TARGET_CLI = $(BIN_DIR)/gssk
TARGET_COMPARE = $(BIN_DIR)/csv_compare
TARGET_SIM = $(BIN_DIR)/giannantoni_sim

# ──────────────────────────────────────────────────────────────
# Containerised Linux toolchains (Apple `container` CLI)
#
# One thing macOS cannot verify locally:
#   gcc  — /usr/bin/gcc here is Apple clang; real GCC emits warnings clang
#          does not, and CFLAGS carries -Werror.
#
# The ubuntu image resolves wrong on Apple silicon, so --platform is
# explicit: without it `container run` fails with "platform linux/arm64"
# even when the image built fine. (WASM needs no container: the WASI SDK
# below runs natively on macOS and Linux.)
# ──────────────────────────────────────────────────────────────
CONTAINER_BIN    := container
CONTAINER_PLATFORM := linux/amd64
IMAGE_LINUX      := gssk-linux
UBUNTU_VERSION   := 24.04
CWORKDIR         := /work
CRUN              = $(CONTAINER_BIN) run --rm --platform $(CONTAINER_PLATFORM) -v $(shell pwd):$(CWORKDIR)

.PHONY: all clean test test-update test-advanced test-price-node test-ratio test-delivered-work test-price-dynamics test-net-energy test-gnp-loop test-node-types test-unknown-keys test-deactivation test-stage-times test-forcing test-wasm test-carrier-api test-edge-flows test-schema check-version demo directories dist \
        shared asan test-asan coverage-build coverage-report coverage-check \
        fuzz-build fuzz-run test-valgrind bench bench-check bench-gen \
        container-start container-image container-image-linux \
        wasm wasi-sdk wasm-toolchain test-wasm-container test-linux test-linux-clang shell-linux ci-local

all: directories $(TARGET_LIB) $(TARGET_CLI) $(TARGET_COMPARE) $(TARGET_SIM)

directories:
	@mkdir -p $(BIN_DIR) $(LIB_DIR) $(DIST_DIR) tests/results tests/expected

# Static Library
$(TARGET_LIB): $(OBJECTS)
	ar rcs $@ $^

# Every object depends on every header.
#
# Without this the pattern rule has no header prerequisites at all, so editing
# a header rebuilds nothing. That is not merely a staleness annoyance: change a
# struct in include/*.h and the objects that were not rebuilt keep the old
# layout, the archive links without complaint, and the binary segfaults on a
# field offset. That is exactly what happened when `gia_node` gained fields --
# engine.o was rebuilt, sim_main.o and validation.o were not.
#
# Coarse on purpose. This project has two headers; per-object dependency
# generation (-MMD -MP) would be more precise and is not worth the machinery.
HEADERS = $(wildcard $(INC_DIR)/*.h)

$(LIB_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# CLI Tool
$(TARGET_CLI): $(SRC_DIR)/main.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) $(TARGET_LIB) -o $@ $(LDFLAGS)

# Test Utility
$(TARGET_COMPARE): $(TEST_DIR)/csv_compare.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

# ──────────────────────────────────────────────────────────────
# Giannantoni Simulation Targets
# ──────────────────────────────────────────────────────────────
#
# The generative engine is a separate binary from bin/gssk, not a subcommand
# of it. The two engines share this build and cJSON; they share no headers.
# src/engine.c deliberately does not include gssk.h -- see the note at the
# top of include/engine.h for why the Relational Space coordinates do not fit
# through an API shaped for a flat double * of scalar storages.
#
# TARGET_SIM itself is defined up with the other TARGET_* variables, because
# `all:` expands its prerequisite list immediately and would otherwise see it
# as empty.

# Default seed graph if MODEL is not passed explicitly
MODEL ?= examples/giannantoni/input.json

# The Giannantoni library units. Defined here, ahead of every rule that uses
# them, because a prerequisite list is expanded when its rule is read.
GIA_UNIT_SRCS = $(SRC_DIR)/engine.c $(SRC_DIR)/validation.c \
                $(SRC_DIR)/projection.c $(SRC_DIR)/idc.c $(SRC_DIR)/mop.c
GIA_OBJS = $(patsubst $(SRC_DIR)/%.c,$(LIB_DIR)/%.o,$(GIA_UNIT_SRCS))

# Simulation objects. sim_main.o carries the entry point, kept out of
# engine.o so tests can link the engine without one.
SIM_OBJS = $(GIA_OBJS) $(LIB_DIR)/sim_main.o

# engine.o, validation.o and sim_main.o are built by the $(LIB_DIR)/%.o
# pattern rule above; they need no rules of their own.

# libgssk.a is linked for cJSON, which lives in it. No GSSK kernel symbol is
# referenced, so the archive contributes cJSON.o and nothing else.
$(TARGET_SIM): $(SIM_OBJS) $(TARGET_LIB)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.PHONY: simulate simulate-linux test-giannantoni demo-giannantoni bench-giannantoni

# Native execution
simulate: directories $(TARGET_SIM)
	@echo "=== Running Giannantoni Simulation (Native) ==="
	@./$(TARGET_SIM) $(MODEL)

# Containerized execution (Apple container CLI)
simulate-linux: container-image-linux
	@echo "=== Running Giannantoni Simulation (Linux Container) ==="
	$(CRUN) $(IMAGE_LINUX) sh -c 'make CC=gcc $(TARGET_SIM) && ./$(TARGET_SIM) $(MODEL)'

# Side-by-side demo of the two modes. Same binary, same flags, no config
# change between the runs -- which is the point. The mode is not something you
# select; it is what the engine turns out to have done, decided afterwards by
# diffing the output graph against the seed.
#
#   input.json       Odum's work gate: sun is the energy, biomass and the
#                    consumer's feedback are controls. Nothing returns to the
#                    sun or the stores, so an emergent quality closes them.
#   closed_loop.json the same gate with a recycler returning to the sun. Every
#                    component that can close already has; the heat sink
#                    never does (ADR 0015), so the run is functional.
demo-giannantoni: directories $(TARGET_SIM)
	@echo
	@echo "### 1. Seed below maximum ordinality -> expect GENERATIVE"
	@./$(TARGET_SIM) examples/giannantoni/input.json --steps 5 --print \
	    --csv $(TEST_DIR)/results/gia_generative.csv \
	    --out $(TEST_DIR)/results/gia_generative.json
	@echo
	@echo "### 2. Seed with nothing left to close -> expect FUNCTIONAL"
	@./$(TARGET_SIM) examples/giannantoni/closed_loop.json --steps 5 --print \
	    --csv $(TEST_DIR)/results/gia_functional.csv \
	    --out $(TEST_DIR)/results/gia_functional.json
	@echo
	@echo "### 3. Feed run 1's own output back in -> expect FUNCTIONAL (fixed point)"
	@./$(TARGET_SIM) $(TEST_DIR)/results/gia_generative.json \
	    --csv $(TEST_DIR)/results/gia_fixpoint.csv \
	    --out $(TEST_DIR)/results/gia_fixpoint.json

# Comparative stress test: the incipient closed form against RK4 and Euler on
# the same ODE, plus the demonstration that the analytic drift psi survives
# dt -> 0 and is therefore not an integration error.
TARGET_BENCH_GIA = $(BIN_DIR)/bench_giannantoni

$(TARGET_BENCH_GIA): bench/bench_giannantoni.c $(LIB_DIR)/engine.o $(TARGET_LIB)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

bench-giannantoni: directories $(TARGET_BENCH_GIA)
	@./$(TARGET_BENCH_GIA)

# Engine unit tests: drift theorem, duet branches, harmony invariants,
# ordinality and the generative step.
TARGET_TEST_GIA = $(BIN_DIR)/test_giannantoni

$(TARGET_TEST_GIA): $(TEST_DIR)/test_giannantoni.c $(GIA_OBJS) $(TARGET_LIB)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

test-giannantoni: directories $(TARGET_TEST_GIA)
	@echo "=== Giannantoni engine tests ==="
	@./$(TARGET_TEST_GIA)

# The Giannantoni kernel's V&V suite: docs/requirements/vv-plan.md §7, every
# test tagged with the requirements it verifies, every oracle an equation or a
# printed number (ADR 0018).
TARGET_TEST_MOP = $(BIN_DIR)/test_mop

$(TARGET_TEST_MOP): $(TEST_DIR)/test_mop.c $(GIA_OBJS) $(TARGET_LIB)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.PHONY: test-mop
test-mop: directories $(TARGET_TEST_MOP)
	@./$(TARGET_TEST_MOP)

# System tests of bin/giannantoni_sim (docs/requirements/vv-plan.md §1): the
# run report, the CSV and the exit status, checked against docs/requirements/icd.md.
.PHONY: test-mop-cli
test-mop-cli: directories $(TARGET_SIM)
	@sh $(TEST_DIR)/mop_cli.sh

# Reentrancy, memory and structure of the Giannantoni units (NFR-REE-001,
# NFR-MEM-001, NFR-ERR-001, NFR-SEP-001; vv-plan.md §7.8).
#
#   test-mop-threads       two models on two threads == each run alone (T-REE-01)
#   test-mop-threads-tsan  the same under ThreadSanitizer; needs a toolchain
#                          with the TSan runtime, so it is local, not in CI
#   test-mop-asan          the Giannantoni test binaries under ASan + LSan +
#                          UBSan, any finding fatal (T-MEM-01); Linux only,
#                          since LeakSanitizer is not available on macOS
#   check-symbols          no writable data, no exit/abort, no gssk.h
#                          (T-REE-02, T-ERR-01, INS-SEP-01)
GIA_SRCS = $(GIA_UNIT_SRCS) $(SRC_DIR)/cJSON.c
TARGET_TEST_THREADS = $(BIN_DIR)/test_mop_threads
SAN_FLAGS = -std=c99 -Iinclude -g -O1 -fno-omit-frame-pointer

$(TARGET_TEST_THREADS): $(TEST_DIR)/test_mop_threads.c $(GIA_OBJS) $(TARGET_LIB)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lpthread

.PHONY: test-mop-threads
test-mop-threads: directories $(TARGET_TEST_THREADS)
	@./$(TARGET_TEST_THREADS)

.PHONY: test-mop-threads-tsan
test-mop-threads-tsan: directories
	$(CC) $(SAN_FLAGS) -fsanitize=thread $(TEST_DIR)/test_mop_threads.c \
	    $(GIA_SRCS) -o $(BIN_DIR)/test_mop_threads_tsan $(LDFLAGS) -lpthread
	@TSAN_OPTIONS=halt_on_error=1 ./$(BIN_DIR)/test_mop_threads_tsan

.PHONY: test-mop-asan
test-mop-asan: directories
	$(CC) $(SAN_FLAGS) -fsanitize=address,undefined -fno-sanitize-recover=all \
	    $(TEST_DIR)/test_giannantoni.c $(GIA_SRCS) \
	    -o $(BIN_DIR)/test_giannantoni_asan $(LDFLAGS)
	$(CC) $(SAN_FLAGS) -fsanitize=address,undefined -fno-sanitize-recover=all \
	    $(TEST_DIR)/test_mop_threads.c $(GIA_SRCS) \
	    -o $(BIN_DIR)/test_mop_threads_asan $(LDFLAGS) -lpthread
	@ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ./$(BIN_DIR)/test_giannantoni_asan > /dev/null 2>&1 \
	    || { echo "test-mop-asan: test_giannantoni FAILED under ASan/LSan/UBSan"; \
	         ASAN_OPTIONS=detect_leaks=1 ./$(BIN_DIR)/test_giannantoni_asan 2>&1 | grep -E 'ERROR|SUMMARY|runtime error|FAIL' | head -20; exit 1; }
	@ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ./$(BIN_DIR)/test_mop_threads_asan
	$(CC) $(SAN_FLAGS) -fsanitize=address,undefined -fno-sanitize-recover=all \
	    $(TEST_DIR)/test_mop.c $(GIA_SRCS) -o $(BIN_DIR)/test_mop_asan $(LDFLAGS)
	@ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 ./$(BIN_DIR)/test_mop_asan > /dev/null 2>&1 \
	    || { echo "test-mop-asan: test_mop FAILED under ASan/LSan/UBSan"; \
	         ASAN_OPTIONS=detect_leaks=1 ./$(BIN_DIR)/test_mop_asan 2>&1 | grep -E 'ERROR|SUMMARY|runtime error|FAIL' | head -20; exit 1; }
	@echo "test-mop-asan: no ASan, LSan or UBSan finding"

.PHONY: check-symbols
check-symbols: all
	@sh scripts/check_symbols.sh

# NFR-API-001 (ADR 0018 rule 3): every function idc.h, relational.h and mop.h
# declare (and IF-API-005's engine.h additions) is called by a test, and no
# unit carries a not-implemented stub. The self-test runs the check against
# icd.md's own API blocks.
.PHONY: check-api-called test-api-called
check-api-called:
	@sh scripts/check_api_called.sh
test-api-called:
	@sh $(TEST_DIR)/api_called_selftest.sh

# Requirements traceability (docs/requirements/README.md): every requirement
# is verified by a catalogued entry, every `implemented` claim is backed by a
# `Verifies:` tag in tests/ or scripts/, and no tag names something that does
# not exist. POSIX sh + awk, so it needs no build and runs anywhere.
.PHONY: check-trace
check-trace:
	@sh scripts/check_trace.sh

# Tests
#
# This glob is deliberately shallow. The Giannantoni engine speaks its own
# model vocabulary -- `system_name`, nodes keyed by `type`, edges by
# `flow_type` -- which the GSSK kernel's schema rejects on sight ("unknown
# top-level key 'system_name'"). That is correct: they are two languages for
# two engines, not one language one of them gets wrong. The generative seeds
# therefore live one level down in examples/giannantoni/, which keeps them
# out of both this glob and the independent one in
# scripts/validate_models.py -- a structural separation rather than an
# exclusion list in each place that would have to be kept in step.
MODELS = $(wildcard examples/*.json)
RESULTS = $(patsubst examples/%.json,tests/results/%.csv,$(MODELS))

# An `X_annotated.json` is documentation, not a second model: it is `X.json`
# with `_`-prefixed commentary the kernel ignores. Its golden CSV therefore has
# to be byte-identical to the plain one, and `make test` checks the two
# trajectories against each other so an edit to one file cannot silently make
# the annotated twin describe a model that is no longer running.
ANNOTATED = $(wildcard examples/*_annotated.json)

# A model with no golden file used to print SKIPPED and pass (PLAN.md §1 B1),
# so a model could "pass" by having no expected output. Now it fails, unless
# tests/skip_allowlist.txt names it with a reason (`name  # why`). An entry for
# a model that does have a golden file is stale and also fails, so the list
# cannot outgrow what it needs. ADR 0018 rule 1; self-tested by
# tests/guard_no_skip.sh. EXPECTED_DIR and SKIP_ALLOWLIST exist for that
# self-test; nothing else should set them.
EXPECTED_DIR   ?= tests/expected
SKIP_ALLOWLIST ?= tests/skip_allowlist.txt

test: all check-version test-schema
	@echo "Running Regression Tests..."
	@mkdir -p tests/results
	@for model in $(MODELS); do \
		name=$$(basename $$model .json); \
		why=$$(awk -v n="$$name" '$$1 == n { r = $$0; sub(/^[^#]*#?[ \t]*/, "", r); print (r == "" ? "-" : r); exit }' $(SKIP_ALLOWLIST) 2>/dev/null); \
		printf "Testing %s... " "$$name"; \
		./bin/gssk $$model tests/results/$$name.csv > /dev/null 2>&1; \
		if [ -f $(EXPECTED_DIR)/$$name.csv ]; then \
			if [ -n "$$why" ]; then echo "FAILED (stale entry: $(SKIP_ALLOWLIST) allowlists $$name, which has a golden file)"; exit 1; fi; \
			./bin/csv_compare $(EXPECTED_DIR)/$$name.csv tests/results/$$name.csv; \
			if [ $$? -eq 0 ]; then echo "PASSED"; \
			else echo "FAILED"; exit 1; fi; \
		elif [ -n "$$why" ] && [ "$$why" != "-" ]; then \
			echo "SKIPPED (allowlisted: $$why)"; \
		elif [ "$$why" = "-" ]; then \
			echo "FAILED (no expected output, and the $(SKIP_ALLOWLIST) entry gives no reason after '#')"; exit 1; \
		else \
			echo "FAILED (no expected output, and $$name is not in $(SKIP_ALLOWLIST); a new kernel model needs 'make test-update', ADR 0018)"; exit 1; \
		fi; \
	done
	@echo "Checking annotated twins against their plain models..."
	@for ann in $(ANNOTATED); do \
		name=$$(basename $$ann .json); \
		plain=$${name%_annotated}; \
		printf "Twin %s == %s... " "$$name" "$$plain"; \
		./bin/csv_compare tests/results/$$plain.csv tests/results/$$name.csv; \
		if [ $$? -eq 0 ]; then echo "PASSED"; \
		else echo "FAILED (an annotated variant has drifted from the model it documents)"; exit 1; fi; \
	done

# The guard above, tested: a removed golden file must turn `make test` red.
.PHONY: test-guard-no-skip
test-guard-no-skip: all
	@sh $(TEST_DIR)/guard_no_skip.sh

# Quick demo — run two models and print the head of each CSV. Pure C: the
# Python plotting step went with the python/ tree, which this fork does not
# carry, and with it the reason this ran in a container.
demo: all
	@echo "=== Decay model (exponential decay, RK4) ==="
	@$(TARGET_CLI) examples/decay_model.json /tmp/gssk_demo_decay.csv
	@head -6 /tmp/gssk_demo_decay.csv
	@echo ""
	@echo "=== Household model (4-carrier ecological-economy) ==="
	@$(TARGET_CLI) examples/household_model.json /tmp/gssk_demo_household.csv
	@head -3 /tmp/gssk_demo_household.csv
	@echo "... ($$(( $$(wc -l < /tmp/gssk_demo_household.csv) - 1 )) data rows, $$(head -1 /tmp/gssk_demo_household.csv | tr ',' '\n' | wc -l | tr -d ' ') columns)"

test-update: all
	@echo "Updating Expected Test Outputs..."
	@mkdir -p tests/expected
	@for model in $(MODELS); do \
		name=$$(basename $$model .json); \
		echo "Generating expected output for $$name"; \
		./bin/gssk $$model tests/expected/$$name.csv; \
	done

# Advanced API test suite (calibration, ensemble, Phase 7 node types)
TARGET_TEST_ADV = $(BIN_DIR)/test_advanced

$(TARGET_TEST_ADV): $(TEST_DIR)/test_advanced.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-advanced: all $(TARGET_TEST_ADV)
	@echo "Running advanced API tests..."
	@./$(TARGET_TEST_ADV)

# Phase C.0 — price_node reference resolution, round-trip, constant fallback
TARGET_TEST_PRICE = $(BIN_DIR)/test_price_node

$(TARGET_TEST_PRICE): $(TEST_DIR)/test_price_node.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-price-node: all $(TARGET_TEST_PRICE)
	@echo "Running price_node tests..."
	@./$(TARGET_TEST_PRICE)

# Phase C.1 — ratio (division) logic: hand-calculated quotient, epsilon floor,
# RK4 vs IDC agreement.
TARGET_TEST_RATIO = $(BIN_DIR)/test_ratio

$(TARGET_TEST_RATIO): $(TEST_DIR)/test_ratio.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-ratio: all $(TARGET_TEST_RATIO)
	@echo "Running ratio logic tests..."
	@./$(TARGET_TEST_RATIO)

# Phase C.2 — delivered work signal: tracking, and non-perturbation of the
# trade it observes.
TARGET_TEST_DW = $(BIN_DIR)/test_delivered_work

$(TARGET_TEST_DW): $(TEST_DIR)/test_delivered_work.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-delivered-work: all $(TARGET_TEST_DW)
	@echo "Running delivered-work tests..."
	@./$(TARGET_TEST_DW)

# Phase C.3 — price dynamics: relaxation toward M/W, the named ratio numerator,
# and the solver paths that only disagree when the Jacobian column is wrong.
TARGET_TEST_PRICEDYN = $(BIN_DIR)/test_price_dynamics

$(TARGET_TEST_PRICEDYN): $(TEST_DIR)/test_price_dynamics.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-price-dynamics: all $(TARGET_TEST_PRICEDYN)
	@echo "Running price dynamics tests..."
	@./$(TARGET_TEST_PRICEDYN)

# Node type validation — an unrecognised node `type` must be an error, not a
# silent fallback to `storage` (ADR 0004). Covers both call sites: GSSK_Init,
# which has full archetype dispatch, and GSSK_AddNode, which has none.
TARGET_TEST_NODETYPE = $(BIN_DIR)/test_node_type_validation

# Limit / threshold logic constants (GIP-0001 G6). The formula was always
# implemented; where C comes from was never written down. These pin the four
# facts the schema now states, including the invisible one: a control-supplied
# C that decays past GSSK_LIMIT_C_EPSILON closes the pathway without an error.
TARGET_TEST_LIMIT = $(BIN_DIR)/test_limit_logic

$(TARGET_TEST_LIMIT): $(TEST_DIR)/test_limit_logic.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

.PHONY: test-limit-logic
test-limit-logic: all $(TARGET_TEST_LIMIT)
	@echo "Running limit/threshold logic tests..."
	@./$(TARGET_TEST_LIMIT)

$(TARGET_TEST_NODETYPE): $(TEST_DIR)/test_node_type_validation.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-node-types: all $(TARGET_TEST_NODETYPE)
	@echo "Running node type validation tests..."
	@./$(TARGET_TEST_NODETYPE)

# Phase C.4 — inflation emerges from net-energy decline. Asserts the CLAIM
# (boom/bust, rising feedback fraction, falling net-per-gross, rising price at
# a constant money supply) rather than a trajectory; `make test` already does
# the golden CSV comparison.
TARGET_TEST_NETENERGY = $(BIN_DIR)/test_net_energy

$(TARGET_TEST_NETENERGY): $(TEST_DIR)/test_net_energy.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-net-energy: all $(TARGET_TEST_NETENERGY)
	@echo "Running net-energy feedback tests..."
	@./$(TARGET_TEST_NETENERGY)

# Phase D.1 — money as a closed conserved GNP loop (Odum Fig. 3). Asserts that
# money's per-carrier conservation error holds within solver_tolerance while
# energy's does not over the same run, and carries the negative control showing
# C.4's open money path scored a perfect 0.0 over no storage at all.
TARGET_TEST_GNPLOOP = $(BIN_DIR)/test_gnp_loop

$(TARGET_TEST_GNPLOOP): $(TEST_DIR)/test_gnp_loop.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-gnp-loop: all $(TARGET_TEST_GNPLOOP)
	@echo "Running GNP-loop conservation tests..."
	@./$(TARGET_TEST_GNPLOOP)

# Forcing functions — one waveform vocabulary, two attachment points. The
# convergence test is the one that catches forcing sampled once per STEP
# instead of once per STAGE; every other test here passes either way.
TARGET_TEST_FORCING = $(BIN_DIR)/test_forcing

# Reversible (barb-less) pathway, GIP-0001 G3 / ADR 0007. Conservation and
# origin/target symmetry are the sharp assertions: any sign or index error in
# build_flow_matrix's four entries breaks them, and neither is visible in a
# golden CSV of node quantities that was regenerated from the same bug.
TARGET_TEST_REVERSIBLE = $(BIN_DIR)/test_reversible

$(TARGET_TEST_REVERSIBLE): $(TEST_DIR)/test_reversible.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

.PHONY: test-reversible
test-reversible: all $(TARGET_TEST_REVERSIBLE)
	@echo "Running reversible pathway tests..."
	@./$(TARGET_TEST_REVERSIBLE)

# n-ary interaction and the subtracting action, GIP-0001 G1 / ADR 0008. The
# sharp assertions are the closed forms and the bit-identity of the n-ary form
# against a collapsed binary one: a build that multiplied only the first
# control still decays, and a golden CSV regenerated from that build agrees
# with itself.
TARGET_TEST_NARY = $(BIN_DIR)/test_interaction_nary

$(TARGET_TEST_NARY): $(TEST_DIR)/test_interaction_nary.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

.PHONY: test-interaction-nary
test-interaction-nary: all $(TARGET_TEST_NARY)
	@echo "Running n-ary interaction / subtract tests..."
	@./$(TARGET_TEST_NARY)

$(TARGET_TEST_FORCING): $(TEST_DIR)/test_forcing.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-forcing: all $(TARGET_TEST_FORCING)
	@echo "Running forcing function tests..."
	@./$(TARGET_TEST_FORCING)

# WASM: the shipped artefact, through the ES module loader, under Node's own
# test runner (no test framework). Three files in tests/wasm/:
#   loader.test.js         the loader's contract: exports, imports, strings, heap
#   corpus.test.js         every regression model, as `make test` runs natively
#   forcing_parity.test.js sin/exp must not silently differ from native
# The native forcing evaluator writes its answers first, for the parity check.
# Needs `make wasm`. The host is not assumed to have Node: use
# `make test-wasm-container`, or NODE=/path/to/node.
NODE ?= node
NODE_IMAGE := docker.io/library/node:22-slim

TARGET_DUMP_FORCING = $(BIN_DIR)/dump_forcing_native

$(TARGET_DUMP_FORCING): $(TEST_DIR)/dump_forcing_native.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-wasm: all $(TARGET_DUMP_FORCING)
	@mkdir -p tests/results
	@./$(TARGET_DUMP_FORCING) tests/results/forcing_native.json
	@test -f $(DIST_DIR)/gssk.wasm -a -f $(DIST_DIR)/gssk.js || { echo "dist/gssk.wasm or dist/gssk.js missing — run 'make wasm' first"; exit 1; }
	@$(NODE) --test "tests/wasm/*.test.js"

# Same, with Node from a container. The native dump runs on the host first,
# so its JSON is the host's libm — exactly what the parity check compares.
test-wasm-container: all $(TARGET_DUMP_FORCING) container-start
	@mkdir -p tests/results
	@./$(TARGET_DUMP_FORCING) tests/results/forcing_native.json
	@test -f $(DIST_DIR)/gssk.wasm -a -f $(DIST_DIR)/gssk.js || { echo "dist/gssk.wasm or dist/gssk.js missing — run 'make wasm' first"; exit 1; }
	$(CONTAINER_BIN) run --rm -v $(shell pwd):$(CWORKDIR) -w $(CWORKDIR) $(NODE_IMAGE) node --test "tests/wasm/*.test.js"

# The "Documentation builds" CI job, run locally with Node from a container.
# It mirrors deploy.yml step for step: Node 20 (the job's node-version, not
# NODE_IMAGE's 22), `npm install`, the VitePress build, then the two checks
# that the theme's design tokens were inlined rather than left as an @import.
# node_modules/ ends up holding Linux builds; it is .gitignored and the host
# never uses it. The token checks run on the host, against the built CSS.
DOCS_NODE_IMAGE := docker.io/library/node:20-slim

.PHONY: docs-build-container
docs-build-container: container-start
	$(CONTAINER_BIN) run --rm -v $(shell pwd):$(CWORKDIR) -w $(CWORKDIR) $(DOCS_NODE_IMAGE) \
		sh -c 'npm install --no-audit --no-fund && npm run docs:build'
	@if grep -rq "@import" docs/.vitepress/dist/assets/*.css; then \
		echo "FAIL: an @import survived the build; the token file did not inline"; exit 1; fi
	@if ! grep -rq -- "--e-accent" docs/.vitepress/dist/assets/*.css; then \
		echo "FAIL: energese tokens are absent from the built docs CSS"; exit 1; fi
	@echo "Docs build OK; design tokens are inlined in the built docs CSS."

# Stage times — the solver must hand each derivative evaluation the right time.
# Pinned BEFORE anything consumes t, so the rest of the suite can hold "nothing
# changed at all" as its criterion. Builds the sources directly with the probe
# macro rather than linking $(TARGET_LIB): the recorder must not exist in the
# shipped library, and this task adds no public API.
TARGET_TEST_STAGET = $(BIN_DIR)/test_stage_times
STAGET_SRC = $(SRC_DIR)/gssk.c $(SRC_DIR)/advanced.c $(SRC_DIR)/cJSON.c

$(TARGET_TEST_STAGET): $(TEST_DIR)/test_stage_times.c $(STAGET_SRC) directories
	$(CC) $(CFLAGS) -DGSSK_STAGE_TIME_PROBE $< $(STAGET_SRC) -o $@ $(LDFLAGS)

test-stage-times: $(TARGET_TEST_STAGET)
	@echo "Running stage time tests..."
	@./$(TARGET_TEST_STAGET)

# Unknown key rejection — a key the parser does not recognise must be an error,
# not silence. Same hazard class as the node-type fallback: a model authored
# against a kernel with a feature this one lacks otherwise runs to completion
# and reports success while producing a different trajectory from the one its
# JSON describes.
TARGET_TEST_UNKNOWNKEY = $(BIN_DIR)/test_unknown_keys

$(TARGET_TEST_UNKNOWNKEY): $(TEST_DIR)/test_unknown_keys.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-unknown-keys: all $(TARGET_TEST_UNKNOWNKEY)
	@echo "Running unknown key tests..."
	@./$(TARGET_TEST_UNKNOWNKEY)

# Deactivation must survive serialise -> reload. GSSK_DeactivateEdge cleared
# `active` AND set k to 0, so the round-trip reproduced the trajectory (k = 0
# kills the flow either way) while losing the flag — and the flag is what the
# ~20 sites that COUNT active elements read: motif detection, the isolated-duet
# test, the closed-system conservation check. GSSK_DeactivateNode was lost
# outright, having no k to hide behind.
TARGET_TEST_DEACT = $(BIN_DIR)/test_deactivation_round_trip

# GSSK_NodeType (GIP-0001 G7). The enum and GSSK_GetNodeTypeString read the
# same field, so they cannot disagree by accident — what these catch is the
# set changing under one and not the other, an ordinal being renumbered (a
# silent break for every WASM consumer), and a composite leaking its own name
# where a primitive belongs.
TARGET_TEST_NODEENUM = $(BIN_DIR)/test_node_type_enum

$(TARGET_TEST_NODEENUM): $(TEST_DIR)/test_node_type_enum.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

.PHONY: test-node-type-enum
test-node-type-enum: all $(TARGET_TEST_NODEENUM)
	@echo "Running node type enum tests..."
	@./$(TARGET_TEST_NODEENUM)

$(TARGET_TEST_DEACT): $(TEST_DIR)/test_deactivation_round_trip.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-deactivation: all $(TARGET_TEST_DEACT)
	@echo "Running deactivation round-trip tests..."
	@./$(TARGET_TEST_DEACT)

# Carrier accessors — the flat getters that keep GSSK_Carrier's struct layout
# from crossing the WASM boundary. Also asserts the flat path and
# GSSK_GetCarrier cannot drift apart.
TARGET_TEST_CARRIER = $(BIN_DIR)/test_carrier_api

$(TARGET_TEST_CARRIER): $(TEST_DIR)/test_carrier_api.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-carrier-api: all $(TARGET_TEST_CARRIER)
	@echo "Running carrier accessor tests..."
	@./$(TARGET_TEST_CARRIER)

# Per-edge flow accessors (GIP-0001 G4). Flow used to be step-local, so a
# consumer could read every node quantity and not one rate. Asserts the cache
# is written at the post-step state, on both step paths, with and without
# quality accounting, and that it agrees with the quality pass it now shares a
# flow expression with.
TARGET_TEST_FLOWS = $(BIN_DIR)/test_edge_flows

$(TARGET_TEST_FLOWS): $(TEST_DIR)/test_edge_flows.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

test-edge-flows: all $(TARGET_TEST_FLOWS)
	@echo "Running per-edge flow accessor tests..."
	@./$(TARGET_TEST_FLOWS)

# Schema conformance — examples/ must match gssk.schema.json.
#
# The kernel does not validate against the schema at load time (ADR 0004), so
# this is what stops the schema and the parser drifting apart. It skips when
# jsonschema is absent rather than failing, so a bare checkout still builds;
# CI installs the dependency so the gate is real there.
# Serialised output is checked too: the schema also has to describe what
# GSSK_SerializeModel/Snapshot emit, which is the format the archival story
# in Phase G rests on. dump_serialized writes those into tests/results.
SER_DIR = tests/results/serialized
TARGET_DUMP_SER = $(BIN_DIR)/dump_serialized

$(TARGET_DUMP_SER): $(TEST_DIR)/dump_serialized.c $(TARGET_LIB)
	$(CC) $(CFLAGS) $< $(TARGET_LIB) -o $@ $(LDFLAGS)

# The npm package version and GSK_VERSION_STRING drifted five majors apart
# before anything noticed, because nothing compared them. Stdlib-only, so it
# runs in any checkout without pip.
check-version:
	@python3 scripts/check_version_sync.py

test-schema: directories $(TARGET_DUMP_SER)
	@rm -rf $(SER_DIR) && mkdir -p $(SER_DIR)
	@./$(TARGET_DUMP_SER) $(SER_DIR) $(MODELS) $(wildcard tests/schema_fixtures/*.json)
	@python3 scripts/validate_models.py

clean:
	rm -rf $(BIN_DIR) $(LIB_DIR) $(DIST_DIR) tests/results coverage/

# ──────────────────────────────────────────────────────────────
# Shared library (for FFI consumers)
# ──────────────────────────────────────────────────────────────
TARGET_SO = $(LIB_DIR)/libgssk.so

shared: directories $(TARGET_SO)

$(TARGET_SO): $(SOURCES)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)

# ──────────────────────────────────────────────────────────────
# AddressSanitizer + UndefinedBehaviorSanitizer build
# Requires clang (CC=clang make asan)
# ──────────────────────────────────────────────────────────────
ASAN_FLAGS = -Wall -Wextra -std=c99 -Iinclude -fPIC -g -O1 \
             -fsanitize=address,undefined -fno-omit-frame-pointer

TARGET_ASAN_CLI = $(BIN_DIR)/gssk_asan

asan: directories
	$(CC) $(ASAN_FLAGS) $(SOURCES) $(SRC_DIR)/main.c -o $(TARGET_ASAN_CLI) $(LDFLAGS)

test-asan: asan $(TARGET_COMPARE)
	@echo "Running regression tests under ASan/UBSan..."
	@mkdir -p tests/results
	@for model in $(MODELS); do \
		name=$$(basename $$model .json); \
		echo -n "ASan $$name... "; \
		if [ ! -f tests/expected/$$name.csv ]; then echo "SKIPPED"; continue; fi; \
		ASAN_OPTIONS=detect_leaks=1 \
		./$(TARGET_ASAN_CLI) $$model tests/results/$$name.csv > /dev/null 2>&1; \
		./bin/csv_compare tests/expected/$$name.csv tests/results/$$name.csv; \
		if [ $$? -eq 0 ]; then echo "PASSED"; else echo "FAILED"; exit 1; fi; \
	done

# ──────────────────────────────────────────────────────────────
# gcov / lcov coverage
# ──────────────────────────────────────────────────────────────
COVERAGE_FLAGS = -Wall -Wextra -std=c99 -Iinclude -fPIC -g -O0 \
                 --coverage -fprofile-arcs -ftest-coverage

TARGET_COV_CLI  = $(BIN_DIR)/gssk_cov
TARGET_COV_ADV  = $(BIN_DIR)/test_advanced_cov
COVERAGE_MIN_LINE = 35

coverage-build: directories
	@mkdir -p coverage
	gcc $(COVERAGE_FLAGS) $(SOURCES) $(SRC_DIR)/main.c -o $(TARGET_COV_CLI) -lm
	gcc $(COVERAGE_FLAGS) $(SOURCES) $(TEST_DIR)/test_advanced.c -o $(TARGET_COV_ADV) -lm

coverage-report: coverage-build
	@mkdir -p coverage/html
	@echo "Running tests to collect coverage data..."
	@mkdir -p tests/results
	@for model in $(MODELS); do \
		name=$$(basename $$model .json); \
		./$(TARGET_COV_CLI) $$model tests/results/$$name.csv > /dev/null 2>&1 || true; \
	done
	@echo "Running advanced API tests for coverage..."
	./$(TARGET_COV_ADV) > /dev/null 2>&1 || true
	lcov --capture --directory . --output-file coverage/lcov.info \
	     --ignore-errors mismatch,unused
	lcov --remove coverage/lcov.info '*/cJSON.c' '*/tests/*' \
	     --output-file coverage/lcov.info \
	     --ignore-errors mismatch,unused
	genhtml coverage/lcov.info --output-directory coverage/html --quiet

# The kernel gate. It used to parse with `grep -oP` (absent from macOS grep)
# and print OK whenever it parsed nothing (PLAN.md §1 B3); an unreadable
# summary now fails. The Giannantoni gate is coverage-gia, below.
coverage-check: coverage-report coverage-gia
	@line_pct=$$(lcov --summary coverage/lcov.info 2>&1 | awk '/lines/ { for (i = 1; i <= NF; i++) if ($$i ~ /^[0-9]+(\.[0-9]+)?%$$/) { sub(/%/, "", $$i); print $$i; exit } }'); \
	if [ -z "$$line_pct" ]; then echo "FAIL: cannot read a line percentage from lcov --summary"; exit 1; fi; \
	echo "Kernel line coverage: $${line_pct}% (gate: $(COVERAGE_MIN_LINE)%)"; \
	if awk -v p="$$line_pct" -v m="$(COVERAGE_MIN_LINE)" 'BEGIN { exit !(p < m) }'; then \
		echo "FAIL: line coverage below $(COVERAGE_MIN_LINE)%"; exit 1; \
	else echo "OK"; fi

# Giannantoni units (NFR-COV-001, ADR 0018 rule 2): >= 90% of lines, measured
# with gcov over every Giannantoni suite. The units are those srs.md names;
# validation.c and projection.c are reported but not gated. A unit listed here
# whose source does not exist yet is skipped, so W2-W8 join the gate by landing.
GIA_COV_MIN    = 90
GIA_COV_UNITS  = engine idc mop relational harmony
GIA_COV_EXTRA  = validation projection
GIA_COV_DIR    = coverage/gia
GCOV          ?= gcov
GIA_COV_FLAGS  = -std=c99 -Iinclude -O0 -g --coverage

.PHONY: coverage-gia
coverage-gia: directories
	@rm -rf $(GIA_COV_DIR) && mkdir -p $(GIA_COV_DIR)
	@for u in $(GIA_COV_UNITS) $(GIA_COV_EXTRA); do \
		[ -f $(SRC_DIR)/$$u.c ] || continue; \
		gcc $(GIA_COV_FLAGS) -c $(SRC_DIR)/$$u.c -o $(GIA_COV_DIR)/$$u.o || exit 1; \
	done
	@gcc -std=c99 -Iinclude -O0 -c $(SRC_DIR)/cJSON.c -o $(GIA_COV_DIR)/cJSON.o
	@objs=$$(ls $(GIA_COV_DIR)/*.o); \
	for t in $(GIA_COV_TESTS); do \
		gcc $(GIA_COV_FLAGS) $(TEST_DIR)/$$t.c $$objs -o $(GIA_COV_DIR)/$$t $(LDFLAGS) -lpthread || exit 1; \
		./$(GIA_COV_DIR)/$$t > $(GIA_COV_DIR)/$$t.log 2>&1 || { echo "coverage-gia: $$t failed"; tail -20 $(GIA_COV_DIR)/$$t.log; exit 1; }; \
	done; \
	gcc $(GIA_COV_FLAGS) $(SRC_DIR)/sim_main.c $$objs -o $(GIA_COV_DIR)/giannantoni_sim $(LDFLAGS) || exit 1; \
	SIM=$(GIA_COV_DIR)/giannantoni_sim sh $(TEST_DIR)/mop_cli.sh > $(GIA_COV_DIR)/mop_cli.log 2>&1 \
		|| { echo "coverage-gia: mop_cli.sh failed"; tail -20 $(GIA_COV_DIR)/mop_cli.log; exit 1; }
	@: > $(GIA_COV_DIR)/gated.txt; : > $(GIA_COV_DIR)/extra.txt
	@for u in $(GIA_COV_UNITS); do \
		[ -f $(SRC_DIR)/$$u.c ] || continue; \
		$(GCOV) -n -o $(GIA_COV_DIR) $(SRC_DIR)/$$u.c 2>/dev/null | awk -v f="$(SRC_DIR)/$$u.c" 'index($$0, "File \047" f "\047") == 1 { p = 1; print; next } p { print; exit }' >> $(GIA_COV_DIR)/gated.txt; \
	done
	@for u in $(GIA_COV_EXTRA); do \
		$(GCOV) -n -o $(GIA_COV_DIR) $(SRC_DIR)/$$u.c 2>/dev/null | awk -v f="$(SRC_DIR)/$$u.c" 'index($$0, "File \047" f "\047") == 1 { p = 1; print; next } p { print; exit }' >> $(GIA_COV_DIR)/extra.txt; \
	done
	@echo "Giannantoni units outside the gate (reported only):"
	@awk '/^File /{ f = $$2 } /^Lines/{ print "  " f " " $$0 }' $(GIA_COV_DIR)/extra.txt
	@echo "Giannantoni units, gated:"
	@sh scripts/coverage_gate.sh $(GIA_COV_MIN) $(GIA_COV_DIR)/gated.txt

# Every Giannantoni test binary built from tests/<name>.c. A W2-W8 suite joins
# coverage by being listed here.
GIA_COV_TESTS = test_giannantoni test_mop_threads test_mop

# The gate's own self-test: garbage, an empty report and 89% must all fail.
.PHONY: test-coverage-gate
test-coverage-gate:
	@sh $(TEST_DIR)/coverage_gate_selftest.sh

# ──────────────────────────────────────────────────────────────
# Valgrind memory-error check (Linux only)
# ──────────────────────────────────────────────────────────────
VALGRIND = valgrind --error-exitcode=1 --leak-check=full \
           --show-leak-kinds=all --track-origins=yes -q

test-valgrind: all
	@echo "Running regression tests under Valgrind..."
	@mkdir -p tests/results
	@for model in $(MODELS); do \
		name=$$(basename $$model .json); \
		echo -n "Valgrind $$name... "; \
		if [ ! -f tests/expected/$$name.csv ]; then echo "SKIPPED"; continue; fi; \
		$(VALGRIND) ./bin/gssk $$model tests/results/$$name.csv > /dev/null 2>&1; \
		if [ $$? -eq 0 ]; then echo "PASSED"; else echo "FAILED (leaks/errors)"; exit 1; fi; \
	done

# ──────────────────────────────────────────────────────────────
# LibFuzzer target (requires clang with -fsanitize=fuzzer)
# ──────────────────────────────────────────────────────────────
FUZZ_FLAGS  = -Wall -std=c99 -Iinclude -g -O1 \
              -fsanitize=fuzzer,address,undefined
TARGET_FUZZ = $(BIN_DIR)/fuzz_gssk
FUZZ_TIMEOUT ?= 30
FUZZ_CORPUS  = tests/fuzz_corpus

fuzz-build: directories
	clang $(FUZZ_FLAGS) $(SOURCES) $(TEST_DIR)/fuzz_gssk.c -lm -o $(TARGET_FUZZ)

fuzz-run: fuzz-build
	@mkdir -p $(FUZZ_CORPUS)
	./$(TARGET_FUZZ) $(FUZZ_CORPUS) -max_total_time=$(FUZZ_TIMEOUT) \
	    -print_final_stats=1 -jobs=1 2>&1 | tail -20

# ──────────────────────────────────────────────────────────────
# Benchmark suite
# ──────────────────────────────────────────────────────────────
BENCH_BASELINE_MS ?= 500    # wall-clock budget for supply_chain_30 (ms)

bench-gen:
	@echo "Generating benchmark models..."
	@python3 bench/gen_bench_models.py

bench: all bench-gen
	@chmod +x bench/run_bench.sh
	@./bench/run_bench.sh

bench-check: all bench-gen
	@chmod +x bench/run_bench.sh
	@./bench/run_bench.sh --regression $(BENCH_BASELINE_MS)

# Sync dist/ without requiring emscripten (schema + TypeScript declarations)
dist: directories
	cp $(SRC_DIR)/gssk.d.ts $(DIST_DIR)/gssk.d.ts
	cp gssk.schema.json $(DIST_DIR)/gssk.schema.json

# ──────────────────────────────────────────────────────────────
# WebAssembly: clang + wasi-libc from the pinned WASI SDK
#
# gssk.wasm is the kernel compiled for wasm32-wasip1 as a "reactor" (a
# library, no main), and src/gssk.js is the whole of its JavaScript side.
# No Emscripten: its generated glue was replaced by that loader, so what
# ships is a standard .wasm and a short, readable ES module.
#
# Release builds are reproduced bit-for-bit, on any architecture, by Guix
# (spike/guix). This SDK is the everyday toolchain, here and in CI.
# `make wasi-sdk` fetches it into tools/ and verifies its SHA-256.
# ──────────────────────────────────────────────────────────────
WASI_SDK_TAG     := wasi-sdk-34
WASI_SDK_VERSION := 34.0
WASI_SDK_HOST    := $(shell uname -m | sed 's/aarch64/arm64/')-$(shell uname -s | tr A-Z a-z | sed 's/darwin/macos/')
WASI_SDK         ?= tools/wasi-sdk-$(WASI_SDK_VERSION)-$(WASI_SDK_HOST)
WASI_SDK_SHA256_arm64-macos  := 9c59398106b417f8f14913380fdf0097a8cc0ff4af9eb3ce0065a859e88d49e9
WASI_SDK_SHA256_x86_64-macos := 87d27fa8adc68dee59bfbf2e22a6d34ef717c34d6bf1d8af2a56fc929d9ce0eb
WASI_SDK_SHA256_arm64-linux  := f7e243dff54d60bcc576e94d6166b69f410f2500ae4a9ceef34315be10e77971
WASI_SDK_SHA256_x86_64-linux := b761e3a0721dbae9c09a0059e5fdb2bf917d1b4a8a7b430fb3b5aafb0984b2c4
SHA256SUM        := $(shell command -v sha256sum >/dev/null 2>&1 && echo sha256sum || echo "shasum -a 256")

# Every symbol the module exports. tests/wasm/loader.test.js reads this
# list, so the build and the test cannot disagree about it.
WASM_EXPORTS := \
	GSSK_Init GSSK_Step GSSK_Reset GSSK_GetState GSSK_GetStateSize GSSK_GetFlows \
	GSSK_GetFlowCount GSSK_GetTStart GSSK_GetTEnd GSSK_GetDt GSSK_GetCurrentTime \
	GSSK_GetStepCount GSSK_GetNodeID GSSK_FindNodeIdx GSSK_GetEdgeID \
	GSSK_FindEdgeIdx GSSK_GetEdgeCount GSSK_GetEdgeK GSSK_SetEdgeK \
	GSSK_GetTransformationRatio GSSK_GetQualityFlow GSSK_GetEdgeQualityFlow \
	GSSK_GetSolverConfidence GSSK_AddNode GSSK_AddEdge GSSK_DeactivateEdge \
	GSSK_DeactivateNode GSSK_ReclassifyNetwork GSSK_SerializeModel \
	GSSK_SerializeSnapshot GSSK_FreeString GSSK_GetSchemaVersion \
	GSSK_GetModelName GSSK_GetModelDescription GSSK_GetModelKernelVersion \
	GSSK_GetModelHash GSSK_GetVersionString GSSK_GetVersionCode \
	GSSK_GetEdgeErrorEstimate GSSK_GetStepErrorEstimate GSSK_GetEventCount \
	GSSK_GetEventTime GSSK_GetEventEdgeID GSSK_GetEventDirection \
	GSSK_EnsembleForecast GSSK_FreeEnsembleResult GSSK_Calibrate \
	GSSK_GetEnsembleNodeCount GSSK_GetEnsembleStepCount GSSK_GetEnsembleMin \
	GSSK_GetEnsembleMax GSSK_GetEnsembleMean GSSK_GetErrorDescription GSSK_Free \
	malloc free GSSK_StepAdaptive GSSK_GetLastStepSize GSSK_GetNextStepSize \
	GSSK_GetConservationError GSSK_SetDiagHooks GSSK_EnableForwardSensitivity \
	GSSK_DisableForwardSensitivity GSSK_GetSensitivity GSSK_RunAdjoint \
	GSSK_GetTransformitySensitivity GSSK_CalibrateGradient \
	GSSK_CalibrateMonteCarlo GSSK_GetMutationCount GSSK_GetMutationRecord \
	GSSK_SetMutationCause GSSK_ClearMutationLog GSSK_ExportMutationLog \
	GSSK_Replay GSSK_GetNodeForcingKind GSSK_GetEdgeForcingKind \
	GSSK_EvaluateNodeForcing GSSK_EvaluateEdgeForcing GSSK_GetCarrierCount \
	GSSK_GetCarrier GSSK_GetNodeCarrier GSSK_GetCarrierID GSSK_GetCarrierUnit \
	GSSK_GetCarrierConserved GSSK_FindCarrierIdx GSSK_GetEdgeCarrier \
	GSSK_GetCarrierConservationError GSSK_GetNodeTypeString GSSK_GetNodeType \
	GSSK_GetArchetypeCount GSSK_GetArchetypeName GSSK_GetCompositeCount \
	GSSK_GetCompositeID GSSK_GetCompositeArchetype GSSK_GetNodeComposite \
	GSSK_GetNodeRole GSSK_GetCompositeMemberCount GSSK_GetCompositeMemberIndex \
	GSSK_SetSeed GSSK_GetSeed GSSK_NextRandom GSSK_NextRandomUniform

# The toolchain. By default the WASI SDK above; guix/gssk.scm overrides these
# to build the release artefact with Guix's clang and its own
# wasi-libc, so this rule stays the one place the WASM build is defined.
WASM_CC              ?= $(WASI_SDK)/bin/clang
WASM_SYSROOT         ?= $(WASI_SDK)/share/wasi-sysroot
WASM_TOOLCHAIN_FLAGS ?=
WASM_TOOLCHAIN_LIBS  ?=

comma := ,
WASM_CFLAGS = --target=wasm32-wasip1 --sysroot=$(WASM_SYSROOT) \
              -mexec-model=reactor -std=c99 -Wall -Wextra -Werror -O3 -Iinclude \
              -Wl,--strip-debug   # the SDK's libc carries DWARF; keep only function names

# The SDK is fetched only when it is the toolchain in use.
wasm: dist $(filter $(WASI_SDK)/bin/clang,$(WASM_CC))
	$(WASM_CC) $(WASM_CFLAGS) $(WASM_TOOLCHAIN_FLAGS) $(SOURCES) \
		$(addprefix -Wl$(comma)--export=,$(WASM_EXPORTS)) $(WASM_TOOLCHAIN_LIBS) -o $(DIST_DIR)/gssk.wasm
	cp $(SRC_DIR)/gssk.js $(DIST_DIR)/gssk.js

wasi-sdk: $(WASI_SDK)/bin/clang

# Rebuild gssk.wasm with a packed toolchain (gssk-toolchain-x86_64-linux.tar.xz,
# from guix.yml run with pack-toolchain), on any x86_64 Linux, without Guix:
#   mkdir tc && tar xf gssk-toolchain-x86_64-linux.tar.xz -C tc
#   tc/bin/make wasm-toolchain TC=$$PWD/tc && sha256sum dist/gssk.wasm
# The result must match the release's gssk-guix.sha256. See guix/README.md.
wasm-toolchain:
	@test -x "$(TC)/bin/clang" || { echo "set TC to the unpacked toolchain directory"; exit 1; }
	PATH="$(TC)/bin:$$PATH" $(MAKE) wasm WASM_CC=clang WASM_SYSROOT=$(TC) \
		WASM_TOOLCHAIN_FLAGS="-nostdlibinc -isystem $(TC)/include/wasm32-wasip1" \
		WASM_TOOLCHAIN_LIBS="-nodefaultlibs -L$(TC)/lib/wasm32-wasip1 -lc $(TC)/lib/wasip1/libclang_rt.builtins-wasm32.a"

$(WASI_SDK)/bin/clang:
	@test -n "$(WASI_SDK_SHA256_$(WASI_SDK_HOST))" || { echo "no pinned WASI SDK for $(WASI_SDK_HOST)"; exit 1; }
	@mkdir -p tools
	curl -fsSL --retry 3 -o tools/wasi-sdk.tar.gz \
		https://github.com/WebAssembly/wasi-sdk/releases/download/$(WASI_SDK_TAG)/wasi-sdk-$(WASI_SDK_VERSION)-$(WASI_SDK_HOST).tar.gz
	echo "$(WASI_SDK_SHA256_$(WASI_SDK_HOST))  tools/wasi-sdk.tar.gz" | $(SHA256SUM) -c -
	tar xzf tools/wasi-sdk.tar.gz -C tools
	rm tools/wasi-sdk.tar.gz

# ──────────────────────────────────────────────────────────────
# Local development
#
# `make dev` serves dev/ — a page that runs any example model through
# dist/gssk.js in the browser — with Vite. gssk.wasm comes from `make wasm`
# (the WASI SDK, natively, in seconds) and Node from the same container image
# test-wasm-container uses, so the host needs only the `container` CLI and
# nothing waits on Guix.
#
#   make dev         http://localhost:5173/
#
# Guix is for releases (guix.yml). To check the page against the bytes a
# release would ship, the Guix path is still here, inside a long-lived
# container that hosts Guix (guix/Containerfile.guix-host):
#
#   make dev-guix    same page, Node and gssk.wasm from the pinned Guix
#                    (first run fetches Guix packages: about 15 minutes)
#   make wasm-guix   the release gssk.wasm into dist/, without the server
#   make guix-down   stop the container; its Guix store is kept
#   make guix-rm     remove it, store and all
#
# The Guix container is long-lived because the Guix store lives in it:
# removing it means fetching every package again.
# ──────────────────────────────────────────────────────────────
GUIX_IMAGE := gssk-guix
GUIX_BOX   := gssk-guix
DEV_PORT   := 5173
GUIX        = guix time-machine -C guix/channels.scm --

.PHONY: guix-image guix-up guix-down guix-rm wasm-guix dev dev-guix

guix-image: container-start
	$(CONTAINER_BIN) build -f guix/Containerfile.guix-host -t $(GUIX_IMAGE) .

guix-up: container-start
	@if ! $(CONTAINER_BIN) inspect $(GUIX_BOX) >/dev/null 2>&1; then \
		$(MAKE) guix-image && \
		$(CONTAINER_BIN) run -d --name $(GUIX_BOX) -c 4 -m 8G -p $(DEV_PORT):$(DEV_PORT) \
			-v $(shell pwd):$(CWORKDIR) $(GUIX_IMAGE) sleep infinity >/dev/null; \
	else \
		$(CONTAINER_BIN) start $(GUIX_BOX) >/dev/null 2>&1 || true; \
	fi

guix-down:
	-$(CONTAINER_BIN) stop $(GUIX_BOX)

guix-rm: guix-down
	-$(CONTAINER_BIN) rm $(GUIX_BOX)

# Store files are read-only; the copies in dist/ must not be, or `make wasm`
# could not overwrite them.
wasm-guix: guix-up dist
	$(CONTAINER_BIN) exec $(GUIX_BOX) sh -c 'cd $(CWORKDIR) && \
		P=$$($(GUIX) build --fallback -f guix/gssk.scm | tail -1) && \
		cp $$P/share/gssk/gssk.wasm $$P/share/gssk/gssk.js $$P/share/gssk/gssk.d.ts $(DIST_DIR)/ && \
		chmod u+w $(DIST_DIR)/gssk.wasm $(DIST_DIR)/gssk.js $(DIST_DIR)/gssk.d.ts && \
		sha256sum $(DIST_DIR)/gssk.wasm'

# npm install runs in the container, so node_modules/ holds Linux builds of
# Vite's native dependencies; it is .gitignored and never used on the host.
VITE = npm install --no-audit --no-fund && \
	npx vite --config dev/vite.config.js --host 0.0.0.0 --port $(DEV_PORT) --strictPort

# -it only from a terminal: without one, `container run -it` refuses to start.
DEV_TTY := $(shell [ -t 0 ] && echo -it)

dev: wasm container-start
	$(CONTAINER_BIN) run --rm $(DEV_TTY) -p $(DEV_PORT):$(DEV_PORT) -v $(shell pwd):$(CWORKDIR) -w $(CWORKDIR) \
		$(NODE_IMAGE) sh -c '$(VITE)'

dev-guix: wasm-guix
	$(CONTAINER_BIN) exec -it $(GUIX_BOX) sh -c 'cd $(CWORKDIR) && \
		$(GUIX) shell -m guix/dev.scm -- sh -c "$(VITE)"'

# ──────────────────────────────────────────────────────────────
# Containerised Linux builds
#
# `make wasm` and a real-GCC build cannot run on macOS directly. These
# targets run them in Linux containers via the Apple `container` CLI.
# Artefacts land in the bind-mounted working tree exactly as a native
# build would.
# ──────────────────────────────────────────────────────────────

# Start the container system daemon (idempotent)
container-start:
	@$(CONTAINER_BIN) system start >/dev/null 2>&1 || true

# Build the native Linux image (matches CI's ubuntu-latest)
container-image-linux: container-start
	$(CONTAINER_BIN) build -f Containerfile.linux -t $(IMAGE_LINUX) \
		--platform $(CONTAINER_PLATFORM) \
		--build-arg UBUNTU_VERSION=$(UBUNTU_VERSION) .

# The build-toolchain image
container-image: container-image-linux

# The suites CI's build-native job runs, plus test-advanced, which the
# pre-push checks require and deploy.yml does not run. One list for both
# compilers: the two recipes used to spell theirs out separately, and drifted --
# gcc missed seven suites CI had gained, clang ran two. Add a suite here when
# you add its step to deploy.yml.
CI_TESTS = test test-advanced test-node-types test-limit-logic test-forcing \
           test-reversible test-interaction-nary test-stage-times \
           test-unknown-keys test-deactivation test-node-type-enum \
           test-carrier-api test-edge-flows test-price-node test-ratio \
           test-delivered-work test-price-dynamics test-net-energy \
           test-gnp-loop test-giannantoni test-mop test-mop-cli test-mop-threads check-symbols \
           test-guard-no-skip test-coverage-gate check-api-called \
           test-api-called check-trace check-version test-schema

# Full native build + CI's suites under real GCC with -Werror.
test-linux: container-image-linux
	$(CRUN) $(IMAGE_LINUX) sh -c 'make clean && make CC=gcc all && make CC=gcc $(CI_TESTS)'

# Same under Linux clang, the other half of CI's build-native matrix.
test-linux-clang: container-image-linux
	$(CRUN) $(IMAGE_LINUX) sh -c 'make clean && make CC=clang all && make CC=clang $(CI_TESTS)'

# Everything CI would catch that macOS cannot: both Linux compilers.
# Leaves the tree holding Linux objects — run `make clean && make all` after.
ci-local: test-linux test-linux-clang
	@echo "──────────────────────────────────────────────"
	@echo "Linux gcc + clang both built."
	@echo "Tree now holds Linux artefacts; run 'make clean && make all' to restore native."

# Interactive shells for debugging a container build
shell-linux: container-image-linux
	$(CONTAINER_BIN) run --rm -it --platform $(CONTAINER_PLATFORM) -v $(shell pwd):$(CWORKDIR) $(IMAGE_LINUX) bash
