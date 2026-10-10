/* Reentrancy of the Giannantoni units (NFR-REE-001, T-REE-01).
 *
 * Two different models are taken through load -> generate -> print on two
 * threads at once, many times over, and every result must equal the one the
 * same model gives when run alone. The oracle is the sequential run of the
 * same code, which is legitimate here because the property under test is
 * "concurrency changes nothing", not any number.
 *
 * The baseline defect this pins (srs NFR-REE-001): gia_generate sorted its
 * open components through a file-scope `static const gia_model *sort_model`
 * read by the qsort comparator, so a second thread could swap the model out
 * from under the first one's sort. Build with -fsanitize=thread to see the
 * race reported directly (`make test-mop-threads-tsan`).
 *
 * `make test-mop-asan` builds this file and tests/test_giannantoni.c with
 * ASan + LSan + UBSan, every finding fatal: that run is T-MEM-01.
 * Verifies: NFR-MEM-001 (T-MEM-01)
 */

#define _POSIX_C_SOURCE 200809L

#include "engine.h"

#include <fcntl.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define ROUNDS 300

/* Below maximum ordinality with several open components, so gia_generate
 * reaches its sort with more than one element. The two models use disjoint
 * ids of different lengths: a comparator reading the other model's node table
 * sorts by the wrong names, or reads past its end. */
static const char *SEED_A =
    "{\"system_name\":\"a\",\"nodes\":["
    " {\"id\":\"src\",\"type\":\"source\",\"initial_value\":2.0},"
    " {\"id\":\"s4\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"s1\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"s3\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"s2\",\"type\":\"storage\",\"current_level\":1.0}],"
    "\"edges\":["
    " {\"source\":\"src\",\"target\":\"s1\",\"weight\":1.0},"
    " {\"source\":\"s1\",\"target\":\"s2\",\"weight\":0.5},"
    " {\"source\":\"s2\",\"target\":\"s3\",\"weight\":0.5},"
    " {\"source\":\"s3\",\"target\":\"s4\",\"weight\":0.5}],"
    "\"simulation_params\":{\"t_val\":1.0,\"derivative_order\":2,"
    " \"generative_mode\":true}}";

static const char *SEED_B =
    "{\"system_name\":\"b\",\"nodes\":["
    " {\"id\":\"zeta_source\",\"type\":\"source\",\"initial_value\":3.0},"
    " {\"id\":\"zeta_hub\",\"type\":\"storage\",\"current_level\":2.0},"
    " {\"id\":\"zeta_c\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"zeta_a\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"zeta_b\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"zeta_d\",\"type\":\"storage\",\"current_level\":1.0},"
    " {\"id\":\"zeta_e\",\"type\":\"storage\",\"current_level\":1.0}],"
    "\"edges\":["
    " {\"source\":\"zeta_source\",\"target\":\"zeta_hub\",\"weight\":1.0},"
    " {\"source\":\"zeta_hub\",\"target\":\"zeta_c\",\"weight\":0.3},"
    " {\"source\":\"zeta_hub\",\"target\":\"zeta_a\",\"weight\":0.3},"
    " {\"source\":\"zeta_hub\",\"target\":\"zeta_b\",\"weight\":0.3},"
    " {\"source\":\"zeta_hub\",\"target\":\"zeta_d\",\"weight\":0.3},"
    " {\"source\":\"zeta_hub\",\"target\":\"zeta_e\",\"weight\":0.3}],"
    "\"simulation_params\":{\"t_val\":1.0,\"derivative_order\":2,"
    " \"generative_mode\":true}}";

/* One complete pass: copy, load, generate, print. Returns a malloc'd string
 * or NULL. The seed is parsed once, before any thread starts, because the
 * vendored cJSON's parser records errors in a process-global (cJSON.c's
 * global_error); that is cJSON's documented limitation, not a Giannantoni
 * unit's, and it is outside NFR-REE-001. */
static char *one_pass(const cJSON *seed) {
    cJSON     *root = cJSON_Duplicate(seed, 1), *out = NULL;
    gia_model  m;
    char      *txt = NULL;

    if (!root) return NULL;
    memset(&m, 0, sizeof(m));
    if (gia_model_load(&m, root)) {
        out = gia_generate(&m);
        if (out) txt = cJSON_PrintUnformatted(out);
        gia_model_free(&m);
    }
    if (out) cJSON_Delete(out);
    cJSON_Delete(root);
    return txt;
}

typedef struct {
    const cJSON *seed;
    const char *expect;
    int         mismatches;
} job;

static void *worker(void *arg) {
    job *j = (job *)arg;
    int  r;
    for (r = 0; r < ROUNDS; r++) {
        char *got = one_pass(j->seed);
        if (!got || strcmp(got, j->expect) != 0) j->mismatches++;
        free(got);
    }
    return NULL;
}

int main(void) {
    char     *ref_a = NULL, *ref_b = NULL;
    cJSON    *seed_a = cJSON_Parse(SEED_A), *seed_b = cJSON_Parse(SEED_B);
    job       ja, jb;
    pthread_t ta, tb;
    int       saved, devnull, failures = 0;

    printf("=== Giannantoni reentrancy (T-REE-01) ===\n");

    /* gia_generate reports its step on stdout; thousands of copies of it are
     * noise here, so stdout is parked on /dev/null while the passes run. */
    fflush(stdout);
    saved   = dup(STDOUT_FILENO);
    devnull = open("/dev/null", O_WRONLY);
    if (saved < 0 || devnull < 0) { perror("dup"); return EXIT_FAILURE; }
    dup2(devnull, STDOUT_FILENO);

    if (seed_a && seed_b) {
        ref_a = one_pass(seed_a);
        ref_b = one_pass(seed_b);
    }

    ja.seed = seed_a; ja.expect = ref_a; ja.mismatches = 0;
    jb.seed = seed_b; jb.expect = ref_b; jb.mismatches = 0;
    if (ref_a && ref_b) {
        pthread_create(&ta, NULL, worker, &ja);
        pthread_create(&tb, NULL, worker, &jb);
        pthread_join(ta, NULL);
        pthread_join(tb, NULL);
    }

    fflush(stdout);
    dup2(saved, STDOUT_FILENO);
    close(saved); close(devnull);

    /* Verifies: NFR-REE-001, BR-011 (T-REE-01) */
    if (!ref_a || !ref_b) {
        printf("  %-58s FAIL\n", "both seeds load and generate sequentially");
        failures++;
    } else {
        /* The seeds must actually exercise the sort: generation adds legs. */
        bool grew_a = strlen(ref_a) > strlen(SEED_A) / 2 &&
                      strstr(ref_a, "max_empower_pathway") != NULL;
        printf("  %-58s %s\n", "seed A is below maximum and grows",
               grew_a ? "PASS" : "FAIL");
        if (!grew_a) failures++;
        printf("  %-58s %s  (%d of %d differ)\n",
               "model A on a thread == model A alone",
               ja.mismatches == 0 ? "PASS" : "FAIL", ja.mismatches, ROUNDS);
        printf("  %-58s %s  (%d of %d differ)\n",
               "model B on a thread == model B alone",
               jb.mismatches == 0 ? "PASS" : "FAIL", jb.mismatches, ROUNDS);
        failures += (ja.mismatches != 0) + (jb.mismatches != 0);
    }
    free(ref_a); free(ref_b);
    cJSON_Delete(seed_a); cJSON_Delete(seed_b);

    printf("\n%s\nfailures: %d\n", failures == 0 ? "ALL PASS" : "FAILURES PRESENT",
           failures);
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
