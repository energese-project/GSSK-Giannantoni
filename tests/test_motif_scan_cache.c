/* K-4: the motif scan is cached, and the cache follows the structure.
 *
 * scan_motifs_internal runs after every step. It enumerates every node pair
 * and triple, O(n^3) with snprintf per triple, and was most of the cost of a
 * step (an 18-node chain spent ~18 us per Euler step in it, ~0.4 us without).
 * Its inputs are the structure alone: the active nodes, their types and the
 * active edges between them. So the kernel now re-enumerates only when those
 * differ from the last full scan, and otherwise keeps the occurrence counts
 * that scan recorded.
 *
 * The risk in a cache is staleness, so every assertion here is about the
 * scan noticing a change. The oracle is the scan's own rules, derived by hand:
 *   - a motif counts toward stability when it occurs >= GSSK_MOTIF_MIN_COUNT
 *     (3) times in a step, and becomes a candidate after
 *     GSSK_MOTIF_MIN_STEPS (10) consecutive such steps;
 *   - the generativity index of a step is new_candidates * mean complexity / dt,
 *     and a 2-node motif with one edge has complexity 1/2.
 * Three parallel source -> storage chains give the source->storage motif an
 * occurrence of exactly 3 per step.
 */

#include "gssk.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static const char *MODEL =
    "{"
    "  \"nodes\": ["
    "    {\"id\": \"s1\", \"type\": \"source\",  \"value\": 10.0},"
    "    {\"id\": \"b1\", \"type\": \"storage\", \"value\": 0.0},"
    "    {\"id\": \"s2\", \"type\": \"source\",  \"value\": 10.0},"
    "    {\"id\": \"b2\", \"type\": \"storage\", \"value\": 0.0},"
    "    {\"id\": \"s3\", \"type\": \"source\",  \"value\": 10.0},"
    "    {\"id\": \"b3\", \"type\": \"storage\", \"value\": 0.0}"
    "  ],"
    "  \"edges\": ["
    "    {\"id\": \"e1\", \"origin\": \"s1\", \"target\": \"b1\", \"logic\": \"linear\", \"params\": {\"k\": 0.1}},"
    "    {\"id\": \"e2\", \"origin\": \"s2\", \"target\": \"b2\", \"logic\": \"linear\", \"params\": {\"k\": 0.1}},"
    "    {\"id\": \"e3\", \"origin\": \"s3\", \"target\": \"b3\", \"logic\": \"linear\", \"params\": {\"k\": 0.1}}"
    "  ],"
    "  \"config\": {\"t_start\": 0, \"t_end\": 100, \"dt\": 1.0, \"method\": \"euler\"}"
    "}";

static int failures = 0;

static void check(int ok, const char *what) {
  printf("  %-64s %s\n", what, ok ? "PASS" : "FAIL");
  if (!ok) failures++;
}

/* The source -> storage motif: the 2-node motif whose types are source and
 * storage. Found by size and types, so the test does not depend on the exact
 * canonical spelling. */
static long find_source_storage(GSSK_Instance *inst) {
  for (size_t i = 0; i < GSSK_GetMotifCount(inst); i++) {
    const char *canon = GSSK_GetMotifCanon(inst, i);
    if (GSSK_GetMotifSize(inst, i) == 2 && canon && strstr(canon, "source") && strstr(canon, "storage"))
      return (long)i;
  }
  return -1;
}

static void step(GSSK_Instance *inst, int times) {
  for (int i = 0; i < times; i++) {
    GSSK_Status st = GSSK_Step(inst, 1.0);
    assert(st == GSSK_SUCCESS || st == GSSK_WARN_SOLVER_DIVERGENCE);
  }
}

static void test_counts_without_change(void) {
  printf("Unchanged structure: the cached counts are the scan's counts\n");
  GSSK_Instance *inst = NULL;
  assert(GSSK_Init(MODEL, &inst) == GSSK_SUCCESS);

  step(inst, 1);
  long m = find_source_storage(inst);
  check(m >= 0, "the source->storage motif is found on the first step");
  check(GSSK_GetMotifOccurrence(inst, (size_t)m) == 3, "it occurs 3 times (three chains)");
  check(GSSK_GetMotifStableSteps(inst, (size_t)m) == 1, "one stable step");

  /* Steps 2..9 reuse the cache: still 3 occurrences, one more stable step each. */
  step(inst, 8);
  check(GSSK_GetMotifOccurrence(inst, (size_t)m) == 3, "still 3 occurrences after 9 steps");
  check(GSSK_GetMotifStableSteps(inst, (size_t)m) == 9, "9 stable steps after 9 steps");
  check(!GSSK_IsMotifCandidate(inst, (size_t)m), "not yet a candidate (needs 10)");

  /* Step 10: becomes a candidate, the only new one: index = 1 * (1/2) / dt. */
  step(inst, 1);
  check(GSSK_IsMotifCandidate(inst, (size_t)m), "a candidate at step 10");
  check(fabs(GSSK_GetGenerativityIndex(inst) - 0.5) < 1e-12, "generativity index 0.5 at step 10");
  step(inst, 1);
  check(GSSK_GetGenerativityIndex(inst) == 0.0, "generativity index 0 once nothing new appears");
  check(GSSK_GetMotifStableSteps(inst, (size_t)m) == 11, "11 stable steps after 11 steps");
  GSSK_Free(inst);
}

static void test_deactivating_an_edge_rescans(void) {
  printf("Deactivating an edge: the next step sees two chains\n");
  GSSK_Instance *inst = NULL;
  assert(GSSK_Init(MODEL, &inst) == GSSK_SUCCESS);
  step(inst, 12);
  long m = find_source_storage(inst);
  assert(m >= 0);

  assert(GSSK_DeactivateEdge(inst, "e3") == GSSK_SUCCESS);
  step(inst, 1);
  check(GSSK_GetMotifOccurrence(inst, (size_t)m) == 2, "occurrence drops to 2");
  check(GSSK_GetMotifStableSteps(inst, (size_t)m) == 0, "below MIN_COUNT: the stable run resets");
  check(GSSK_IsMotifCandidate(inst, (size_t)m), "a candidate stays a candidate");

  /* A third chain again, through a new edge: counted again from the next step. */
  assert(GSSK_AddEdge(inst, "{\"id\":\"e3b\",\"origin\":\"s3\",\"target\":\"b3\",\"logic\":\"linear\",\"params\":{\"k\":0.1}}") ==
         GSSK_SUCCESS);
  step(inst, 1);
  check(GSSK_GetMotifOccurrence(inst, (size_t)m) == 3, "an added edge brings it back to 3");
  check(GSSK_GetMotifStableSteps(inst, (size_t)m) == 1, "and the stable run starts again");
  GSSK_Free(inst);
}

static void test_deactivating_a_node_rescans(void) {
  printf("Deactivating a node: its motifs are no longer counted\n");
  GSSK_Instance *inst = NULL;
  assert(GSSK_Init(MODEL, &inst) == GSSK_SUCCESS);
  step(inst, 3);
  long m = find_source_storage(inst);
  assert(m >= 0);
  assert(GSSK_DeactivateNode(inst, "b1") == GSSK_SUCCESS);
  step(inst, 1);
  check(GSSK_GetMotifOccurrence(inst, (size_t)m) == 2, "occurrence drops to 2 without b1");
  GSSK_Free(inst);
}

static void test_new_structure_adds_motifs(void) {
  printf("A new node and edge: a new motif is recorded\n");
  GSSK_Instance *inst = NULL;
  assert(GSSK_Init(MODEL, &inst) == GSSK_SUCCESS);
  step(inst, 2);
  size_t before = GSSK_GetMotifCount(inst);
  assert(GSSK_AddNode(inst, "{\"id\":\"b4\",\"type\":\"storage\",\"value\":0.0}") == GSSK_SUCCESS);
  assert(GSSK_AddEdge(inst, "{\"id\":\"e4\",\"origin\":\"b1\",\"target\":\"b4\",\"logic\":\"linear\",\"params\":{\"k\":0.1}}") ==
         GSSK_SUCCESS);
  step(inst, 1);
  check(GSSK_GetMotifCount(inst) > before, "the storage->storage motifs appear");
  GSSK_Free(inst);
}

static void test_instances_do_not_share_the_cache(void) {
  printf("Two instances: each scan has its own cache\n");
  GSSK_Instance *a = NULL;
  GSSK_Instance *b = NULL;
  assert(GSSK_Init(MODEL, &a) == GSSK_SUCCESS);
  assert(GSSK_Init(MODEL, &b) == GSSK_SUCCESS);
  step(a, 2);
  step(b, 2);
  assert(GSSK_DeactivateEdge(b, "e1") == GSSK_SUCCESS);
  step(a, 1);
  step(b, 1);
  long ma = find_source_storage(a);
  long mb = find_source_storage(b);
  check(GSSK_GetMotifOccurrence(a, (size_t)ma) == 3, "the untouched instance still counts 3");
  check(GSSK_GetMotifOccurrence(b, (size_t)mb) == 2, "the changed one counts 2");
  GSSK_Free(a);
  GSSK_Free(b);
}

int main(void) {
  test_counts_without_change();
  test_deactivating_an_edge_rescans();
  test_deactivating_a_node_rescans();
  test_new_structure_adds_motifs();
  test_instances_do_not_share_the_cache();
  printf(failures ? "FAILED (%d failures)\n" : "OK (0 failures)\n", failures);
  return failures ? 1 : 0;
}
