/* libFuzzer harness for seams_system_set_atoms, seams_system_set_positions
 * and seams_system_stage_positions.
 *
 * The input picks the atom count (0 to 64), which pointers are NULL,
 * wrong lengths and counts, NaN and infinite values, raw or bounded
 * coordinates, the staging mode, and up to eight reloads. After every
 * call the harness reads the frame back:
 *   success   the stored coordinates equal the input bit for bit, and a
 *             reload keeps the types and ids
 *   failure   a message is set, and coordinates, cell, types and ids are
 *             those stored before the call
 * Any difference aborts.
 */
#include "seams.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define N_MAX 64

struct reader {
  const uint8_t *data;
  size_t size;
  size_t at;
};

static uint8_t byte(struct reader *in) {
  return in->at < in->size ? in->data[in->at++] : 0;
}

static double number(struct reader *in, int raw) {
  if (raw) {
    uint64_t bits = 0;
    for (int k = 0; k < 8; ++k) {
      bits |= (uint64_t)byte(in) << (8 * k);
    }
    double value;
    memcpy(&value, &bits, sizeof value);
    return value;
  }
  const int16_t step = (int16_t)(byte(in) | (byte(in) << 8));
  return step / 64.0;
}

struct snapshot {
  int ok;
  size_t n;
  double xyz[3 * N_MAX];
  double box[9];
  int types[N_MAX];
  int ids[N_MAX];
};

static void take(seams_system *sys, struct snapshot *shot) {
  memset(shot, 0, sizeof *shot);
  if (seams_system_natoms(sys, &shot->n) != SEAMS_OK || shot->n == 0) {
    return;
  }
  if (shot->n > N_MAX) {
    abort();
  }
  shot->ok = seams_system_positions(sys, shot->xyz, 3 * shot->n) == SEAMS_OK &&
             seams_system_box(sys, shot->box, 9) == SEAMS_OK &&
             seams_system_types(sys, shot->types, shot->n) == SEAMS_OK &&
             seams_system_ids(sys, shot->ids, shot->n) == SEAMS_OK;
  if (!shot->ok) {
    abort();
  }
}

static int same(const struct snapshot *a, const struct snapshot *b) {
  return a->ok == b->ok && a->n == b->n &&
         memcmp(a->xyz, b->xyz, sizeof(double) * 3 * a->n) == 0 &&
         memcmp(a->box, b->box, sizeof a->box) == 0 &&
         memcmp(a->types, b->types, sizeof(int) * a->n) == 0 &&
         memcmp(a->ids, b->ids, sizeof(int) * a->n) == 0;
}

static void failed(seams_status status, seams_system *sys, const struct snapshot *before) {
  if (status == SEAMS_OK) {
    return;
  }
  const char *message = seams_last_error_message();
  if (message == NULL || message[0] == '\0') {
    abort();
  }
  (void)seams_handle_error(sys);
  struct snapshot after;
  take(sys, &after);
  if (!same(before, &after)) {
    abort();
  }
}

static void fill(struct reader *in, double *xyz, size_t count, int raw) {
  for (size_t k = 0; k < count; ++k) {
    xyz[k] = number(in, raw);
  }
  const uint8_t poison = byte(in);
  if (count > 0 && (poison & 0xc0) != 0) {
    xyz[poison % count] = (poison & 0x40) ? NAN : ((poison & 0x20) ? -INFINITY : INFINITY);
  }
}

static void cell(struct reader *in, double *box, int raw) {
  for (int k = 0; k < 9; ++k) {
    box[k] = number(in, raw);
  }
  if (!raw) {
    for (int k = 3; k < 6; ++k) {
      box[k] = 4.0 + fabs(box[k]);
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  struct reader in = {data, size, 0};
  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    return 0;
  }
  const uint8_t mode = byte(&in);
  const size_t n = byte(&in) % (N_MAX + 1);
  double xyz[3 * N_MAX + 3];
  double box[9];
  int types[N_MAX];
  int ids[N_MAX];
  fill(&in, xyz, 3 * n, mode & 1);
  cell(&in, box, mode & 2);
  for (size_t i = 0; i < n; ++i) {
    types[i] = 1 + byte(&in) % 3;
    ids[i] = (mode & 4) ? (int)(n - i) : (int)i + 1;
  }
  if (mode & 8) {
    (void)seams_system_stage_positions(sys, 1);
  }

  struct snapshot before, after;
  take(sys, &before);
  seams_status status = seams_system_set_atoms(sys, (mode & 16) ? NULL : xyz, n,
                                               (mode & 32) ? NULL : box, (mode & 64) ? 8 : 9,
                                               (mode & 128) ? NULL : types,
                                               (mode & 4) ? ids : NULL);
  failed(status, sys, &before);
  if (status == SEAMS_OK) {
    take(sys, &after);
    if (after.n != n || memcmp(after.xyz, xyz, sizeof(double) * 3 * n) != 0) {
      abort();
    }
  }

  const int reloads = 1 + byte(&in) % 8;
  for (int r = 0; r < reloads && in.at < in.size; ++r) {
    const uint8_t how = byte(&in);
    const uint8_t extra = byte(&in);
    if (how & 1) {
      (void)seams_system_stage_positions(sys, (how & 2) != 0);
    }
    if (extra & 1) {
      const size_t m = (extra & 2) && n > 0 ? n - 1 : n;
      for (size_t i = 0; i < n; ++i) {
        types[i] = 1 + types[i] % 3;
      }
      take(sys, &before);
      status = seams_system_set_atoms(sys, xyz, m, box, 9, types, NULL);
      failed(status, sys, &before);
      if (status == SEAMS_OK) {
        take(sys, &after);
        if (after.n != m || memcmp(after.xyz, xyz, sizeof(double) * 3 * m) != 0 ||
            memcmp(after.types, types, sizeof(int) * m) != 0) {
          abort();
        }
      }
    }
    take(sys, &before);
    const size_t stored = before.n;
    size_t count = stored;
    if (how & 4) {
      count = (how & 8) ? stored + 1 : (stored > 0 ? stored - 1 : 0);
    }
    double next[3 * N_MAX + 3];
    double next_box[9];
    fill(&in, next, 3 * count, how & 16);
    cell(&in, next_box, how & 32);
    const int with_box = (how & 64) != 0;
    seams_system *target = (how & 128) && (how & 2) ? NULL : sys;
    status = seams_system_set_positions(target, (how & 128) && !(how & 2) ? NULL : next, count,
                                        with_box ? next_box : NULL, with_box ? 9 : 0);
    if (target == NULL) {
      if (status != SEAMS_ERR_NULL) {
        abort();
      }
      take(sys, &after);
      if (!same(&before, &after)) {
        abort();
      }
      continue;
    }
    failed(status, sys, &before);
    if (status == SEAMS_OK) {
      take(sys, &after);
      if (after.n != count || memcmp(after.xyz, next, sizeof(double) * 3 * count) != 0 ||
          memcmp(after.types, before.types, sizeof(int) * count) != 0 ||
          memcmp(after.ids, before.ids, sizeof(int) * count) != 0) {
        abort();
      }
      if (!with_box && memcmp(after.box, before.box, sizeof after.box) != 0) {
        abort();
      }
      if (how & 8) {
        int labels[N_MAX];
        (void)seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                                 count, NULL);
      }
    }
  }
  if (seams_system_stage_positions(NULL, 1) != SEAMS_ERR_NULL) {
    abort();
  }
  seams_system_destroy(sys);
  return 0;
}
