/* Copyright (c) d-SEAMS developers
 * SPDX-License-Identifier: MIT
 *
 * Lua 5.4 and LuaJIT module over include/seams.h.
 * Handles are userdata. __gc and __close release them.
 * A status other than SEAMS_OK is a Lua error whose text is
 * seams_last_error when that string is set, otherwise the status name.
 */

#include "seams.h"

#include <lauxlib.h>
#include <lua.h>
#include <stdlib.h>
#include <string.h>

#if defined(__has_include)
#if __has_include(<luajit.h>)
#include <luajit.h>
#endif
#endif

#if LUA_VERSION_NUM < 502
#define luaL_len(L, i) ((lua_Integer)lua_objlen((L), (i)))
#endif

#if LUA_VERSION_NUM < 502 && !defined(LUAJIT_VERSION)
static void luaL_setfuncs(lua_State *L, const luaL_Reg *reg, int nup) {
  for (; reg->name != NULL; reg++) {
    int i;
    for (i = 0; i < nup; i++) {
      lua_pushvalue(L, -nup);
    }
    lua_pushcclosure(L, reg->func, nup);
    lua_setfield(L, -(nup + 2), reg->name);
  }
  lua_pop(L, nup);
}
#ifndef luaL_newlib
#define luaL_newlib(L, reg) (lua_newtable(L), luaL_setfuncs((L), (reg), 0))
#endif
#endif

typedef struct {
  seams_system *sys;
} SystemBox;

typedef struct {
  seams_neighbours *nb;
} NeighbourBox;

typedef struct {
  seams_rings *rings;
} RingsBox;

typedef struct {
  seams_cage_set *cages;
} CageBox;

static int raise_status(lua_State *L, seams_status status) {
  const char *detail = seams_last_error();
  if (detail != NULL && detail[0] != '\0') {
    return luaL_error(L, "%s", detail);
  }
  detail = seams_status_string(status);
  return luaL_error(L, "%s", detail != NULL ? detail : "seams error");
}

static int check_ok(lua_State *L, seams_status status) {
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  return 0;
}

static SystemBox *system_box(lua_State *L, int index) {
  SystemBox *box = (SystemBox *)luaL_checkudata(L, index, "dseams.System");
  if (box->sys == NULL) {
    luaL_error(L, "system is closed");
  }
  return box;
}

static NeighbourBox *neighbour_box(lua_State *L, int index) {
  NeighbourBox *box = (NeighbourBox *)luaL_checkudata(L, index, "dseams.Neighbours");
  if (box->nb == NULL) {
    luaL_error(L, "neighbours are closed");
  }
  return box;
}

static RingsBox *rings_box(lua_State *L, int index) {
  RingsBox *box = (RingsBox *)luaL_checkudata(L, index, "dseams.Rings");
  if (box->rings == NULL) {
    luaL_error(L, "rings are closed");
  }
  return box;
}

/* 0 on success. On failure pushes an error string and leaves *out freed. */
static int read_numbers(lua_State *L, int index, double **out, int *count) {
  int n;
  int i;
  double *buf;
  *out = NULL;
  *count = 0;
  if (!lua_istable(L, index)) {
    lua_pushstring(L, "expected a table of numbers");
    return 1;
  }
  n = (int)luaL_len(L, index);
  buf = (double *)malloc((size_t)(n > 0 ? n : 1) * sizeof(double));
  if (buf == NULL) {
    lua_pushstring(L, "out of memory");
    return 1;
  }
  for (i = 0; i < n; i++) {
    lua_rawgeti(L, index, i + 1);
    if (!lua_isnumber(L, -1)) {
      lua_pop(L, 1);
      free(buf);
      lua_pushstring(L, "table entry is not a number");
      return 1;
    }
    buf[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  *out = buf;
  *count = n;
  return 0;
}

static int read_ints(lua_State *L, int index, int **out, int expect) {
  int n;
  int i;
  int *buf;
  *out = NULL;
  if (lua_isnoneornil(L, index)) {
    return 0;
  }
  if (!lua_istable(L, index)) {
    lua_pushstring(L, "expected a table of integers");
    return 1;
  }
  n = (int)luaL_len(L, index);
  if (expect >= 0 && n != expect) {
    lua_pushfstring(L, "integer table length is %d, want %d", n, expect);
    return 1;
  }
  buf = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
  if (buf == NULL) {
    lua_pushstring(L, "out of memory");
    return 1;
  }
  for (i = 0; i < n; i++) {
    lua_rawgeti(L, index, i + 1);
    if (!lua_isnumber(L, -1)) {
      lua_pop(L, 1);
      free(buf);
      lua_pushstring(L, "table entry is not a number");
      return 1;
    }
    buf[i] = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
  *out = buf;
  return 0;
}

static int fail_msg(lua_State *L) {
  const char *msg = lua_tostring(L, -1);
  return luaL_error(L, "%s", msg != NULL ? msg : "seams error");
}

static void push_ints(lua_State *L, const int *values, int n) {
  int i;
  lua_createtable(L, n, 0);
  for (i = 0; i < n; i++) {
    lua_pushinteger(L, values[i]);
    lua_rawseti(L, -2, i + 1);
  }
}

static void push_doubles(lua_State *L, const double *values, int n) {
  int i;
  lua_createtable(L, n, 0);
  for (i = 0; i < n; i++) {
    lua_pushnumber(L, values[i]);
    lua_rawseti(L, -2, i + 1);
  }
}

static int sys_close(lua_State *L) {
  SystemBox *box = (SystemBox *)luaL_checkudata(L, 1, "dseams.System");
  seams_system_destroy(box->sys);
  box->sys = NULL;
  return 0;
}

static int nb_close(lua_State *L) {
  NeighbourBox *box = (NeighbourBox *)luaL_checkudata(L, 1, "dseams.Neighbours");
  seams_neighbours_destroy(box->nb);
  box->nb = NULL;
  return 0;
}

static int rings_close(lua_State *L) {
  RingsBox *box = (RingsBox *)luaL_checkudata(L, 1, "dseams.Rings");
  seams_rings_destroy(box->rings);
  box->rings = NULL;
  return 0;
}

static int sys_set_frame(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double *xyz = NULL;
  double *cell = NULL;
  int *types = NULL;
  int *ids = NULL;
  int nxyz = 0;
  int ncell = 0;
  int atoms;
  seams_status status;
  if (read_numbers(L, 2, &xyz, &nxyz) || read_numbers(L, 3, &cell, &ncell)) {
    free(xyz);
    free(cell);
    return fail_msg(L);
  }
  if (nxyz <= 0 || nxyz % 3 != 0) {
    free(xyz);
    free(cell);
    return luaL_error(L, "coordinates must be a flat x,y,z table");
  }
  if (ncell != 9) {
    free(xyz);
    free(cell);
    return luaL_error(L, "box must be nine numbers");
  }
  atoms = nxyz / 3;
  if (read_ints(L, 4, &types, atoms) || read_ints(L, 5, &ids, atoms)) {
    free(xyz);
    free(cell);
    free(types);
    free(ids);
    return fail_msg(L);
  }
  status = seams_system_set_frame(box->sys, xyz, atoms, cell, types, ids);
  free(xyz);
  free(cell);
  free(types);
  free(ids);
  check_ok(L, status);
  lua_pushinteger(L, seams_system_count(box->sys));
  return 1;
}

static int sys_set_ortho(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double *xyz = NULL;
  int *types = NULL;
  int *ids = NULL;
  int nxyz = 0;
  int atoms;
  double ox, oy, oz, lx, ly, lz;
  seams_status status;
  if (read_numbers(L, 2, &xyz, &nxyz)) {
    return fail_msg(L);
  }
  if (nxyz <= 0 || nxyz % 3 != 0) {
    free(xyz);
    return luaL_error(L, "coordinates must be a flat x,y,z table");
  }
  atoms = nxyz / 3;
  ox = luaL_checknumber(L, 3);
  oy = luaL_checknumber(L, 4);
  oz = luaL_checknumber(L, 5);
  lx = luaL_checknumber(L, 6);
  ly = luaL_checknumber(L, 7);
  lz = luaL_checknumber(L, 8);
  if (read_ints(L, 9, &types, atoms) || read_ints(L, 10, &ids, atoms)) {
    free(xyz);
    free(types);
    free(ids);
    return fail_msg(L);
  }
  status = seams_system_set_ortho(box->sys, xyz, atoms, ox, oy, oz, lx, ly, lz, types, ids);
  free(xyz);
  free(types);
  free(ids);
  check_ok(L, status);
  lua_pushinteger(L, seams_system_count(box->sys));
  return 1;
}

static int sys_count(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  lua_pushinteger(L, seams_system_count(box->sys));
  return 1;
}

static int sys_chill_plus(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cutoff = luaL_checknumber(L, 2);
  int type = (int)luaL_checkinteger(L, 3);
  int graph = (int)luaL_checkinteger(L, 4);
  int rule = (int)luaL_checkinteger(L, 5);
  int n = seams_system_count(box->sys);
  int counts[9];
  int *labels;
  seams_status status;
  if (n <= 0) {
    return luaL_error(L, "system has no frame");
  }
  labels = (int *)calloc((size_t)n, sizeof(int));
  if (labels == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_chill_plus(box->sys, cutoff, type, graph, rule, labels, counts);
  if (status != SEAMS_OK) {
    free(labels);
    return raise_status(L, status);
  }
  push_ints(L, labels, n);
  push_ints(L, counts, 9);
  free(labels);
  return 2;
}

static int sys_q4q6(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cutoff = luaL_checknumber(L, 2);
  int type = (int)luaL_checkinteger(L, 3);
  int n = seams_system_count(box->sys);
  double *q4;
  double *q6;
  seams_status status;
  if (n <= 0) {
    return luaL_error(L, "system has no frame");
  }
  q4 = (double *)calloc((size_t)n, sizeof(double));
  q6 = (double *)calloc((size_t)n, sizeof(double));
  if (q4 == NULL || q6 == NULL) {
    free(q4);
    free(q6);
    return luaL_error(L, "out of memory");
  }
  status = seams_q4q6(box->sys, cutoff, type, q4, q6);
  if (status != SEAMS_OK) {
    free(q4);
    free(q6);
    return raise_status(L, status);
  }
  push_doubles(L, q4, n);
  push_doubles(L, q6, n);
  free(q4);
  free(q6);
  return 2;
}

static void push_neighbours(lua_State *L, seams_neighbours *nb) {
  NeighbourBox *box = (NeighbourBox *)lua_newuserdata(L, sizeof(NeighbourBox));
  box->nb = nb;
  luaL_getmetatable(L, "dseams.Neighbours");
  lua_setmetatable(L, -2);
}

static void push_rings(lua_State *L, seams_rings *rings) {
  RingsBox *box = (RingsBox *)lua_newuserdata(L, sizeof(RingsBox));
  box->rings = rings;
  luaL_getmetatable(L, "dseams.Rings");
  lua_setmetatable(L, -2);
}

static int sys_neighbours(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cutoff = luaL_checknumber(L, 2);
  int type = (int)luaL_checkinteger(L, 3);
  int graph = (int)luaL_checkinteger(L, 4);
  seams_neighbours *nb = NULL;
  seams_status status = seams_neighbours_build(box->sys, cutoff, type, graph, &nb);
  check_ok(L, status);
  push_neighbours(L, nb);
  return 1;
}

static int sys_rings(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cutoff = luaL_checknumber(L, 2);
  int type = (int)luaL_checkinteger(L, 3);
  int graph = (int)luaL_checkinteger(L, 4);
  int depth = (int)luaL_checkinteger(L, 5);
  seams_rings *rings = NULL;
  seams_status status = seams_rings_find(box->sys, cutoff, type, graph, depth, &rings);
  check_ok(L, status);
  push_rings(L, rings);
  return 1;
}

static int nb_count(lua_State *L) {
  NeighbourBox *box = neighbour_box(L, 1);
  lua_pushinteger(L, seams_neighbours_count(box->nb));
  return 1;
}

static int nb_csr(lua_State *L) {
  NeighbourBox *box = neighbour_box(L, 1);
  size_t n_offsets = 0;
  size_t nnz = 0;
  int *offsets = NULL;
  int *cols = NULL;
  seams_status status =
      seams_neighbours_copy(box->nb, NULL, 0, &n_offsets, NULL, 0, &nnz);
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  offsets = (int *)calloc(n_offsets > 0 ? n_offsets : 1, sizeof(int));
  cols = (int *)calloc(nnz > 0 ? nnz : 1, sizeof(int));
  if (offsets == NULL || cols == NULL) {
    free(offsets);
    free(cols);
    return luaL_error(L, "out of memory");
  }
  status = seams_neighbours_copy(box->nb, offsets, n_offsets, &n_offsets, cols, nnz, &nnz);
  if (status != SEAMS_OK) {
    free(offsets);
    free(cols);
    return raise_status(L, status);
  }
  push_ints(L, offsets, (int)n_offsets);
  push_ints(L, cols, (int)nnz);
  lua_pushinteger(L, (lua_Integer)nnz);
  free(offsets);
  free(cols);
  return 3;
}

static int rings_count(lua_State *L) {
  RingsBox *box = rings_box(L, 1);
  lua_pushinteger(L, seams_rings_count(box->rings));
  return 1;
}

static int rings_histogram(lua_State *L) {
  RingsBox *box = rings_box(L, 1);
  int len = (int)luaL_checkinteger(L, 2);
  int *counts;
  seams_status status;
  if (len <= 0) {
    return luaL_error(L, "histogram length must be positive");
  }
  counts = (int *)calloc((size_t)len, sizeof(int));
  if (counts == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_rings_histogram(box->rings, counts, len);
  if (status != SEAMS_OK) {
    free(counts);
    return raise_status(L, status);
  }
  push_ints(L, counts, len);
  free(counts);
  return 1;
}

static int rings_labels(lua_State *L) {
  RingsBox *box = rings_box(L, 1);
  int n = seams_rings_count(box->rings);
  int *sizes = NULL;
  int *symmetry = NULL;
  seams_status status;
  if (n > 0) {
    sizes = (int *)calloc((size_t)n, sizeof(int));
    symmetry = (int *)calloc((size_t)n, sizeof(int));
    if (sizes == NULL || symmetry == NULL) {
      free(sizes);
      free(symmetry);
      return luaL_error(L, "out of memory");
    }
  }
  status = seams_rings_labels(box->rings, sizes, symmetry);
  if (status != SEAMS_OK) {
    free(sizes);
    free(symmetry);
    return raise_status(L, status);
  }
  if (n == 0) {
    lua_createtable(L, 0, 0);
    lua_createtable(L, 0, 0);
  } else {
    push_ints(L, sizes, n);
    push_ints(L, symmetry, n);
  }
  free(sizes);
  free(symmetry);
  return 2;
}

static int rings_fraction(lua_State *L) {
  RingsBox *box = rings_box(L, 1);
  int size = (int)luaL_checkinteger(L, 2);
  double fraction = 0.0;
  seams_status status = seams_rings_fraction(box->rings, size, &fraction);
  check_ok(L, status);
  lua_pushnumber(L, fraction);
  return 1;
}

static int l_system(lua_State *L) {
  SystemBox *box = (SystemBox *)lua_newuserdata(L, sizeof(SystemBox));
  box->sys = seams_system_create();
  luaL_getmetatable(L, "dseams.System");
  lua_setmetatable(L, -2);
  if (box->sys == NULL) {
    return raise_status(L, SEAMS_ERR_NOMEM);
  }
  return 1;
}

static int l_version(lua_State *L) {
  lua_pushstring(L, seams_version_string());
  return 1;
}

static int l_abi(lua_State *L) {
  lua_pushinteger(L, seams_abi());
  return 1;
}

static int l_status_string(lua_State *L) {
  lua_pushstring(L, seams_status_string((seams_status)luaL_checkinteger(L, 1)));
  return 1;
}

static int l_label_name(lua_State *L) {
  lua_pushstring(L, seams_label_name((int)luaL_checkinteger(L, 1)));
  return 1;
}

static int l_ring_name(lua_State *L) {
  lua_pushstring(L, seams_ring_symmetry_name((int)luaL_checkinteger(L, 1)));
  return 1;
}

static int l_set_threads(lua_State *L) {
  seams_status status = seams_set_num_threads((int)luaL_checkinteger(L, 1));
  check_ok(L, status);
  lua_pushinteger(L, seams_get_num_threads());
  return 1;
}

static int l_threads(lua_State *L) {
  lua_pushinteger(L, seams_get_num_threads());
  return 1;
}

static int l_last_error(lua_State *L) {
  const char *detail = seams_last_error();
  lua_pushstring(L, detail != NULL ? detail : "");
  return 1;
}

static int sys_read(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  const char *path = luaL_checkstring(L, 2);
  int frame = (int)luaL_checkinteger(L, 3);
  int type_filter = (int)luaL_checkinteger(L, 4);
  check_ok(L, seams_system_read(box->sys, path, frame, type_filter));
  lua_pushinteger(L, seams_system_count(box->sys));
  return 1;
}

static int l_nframes(lua_State *L) {
  size_t nframes = 0;
  check_ok(L, seams_path_nframes(luaL_checkstring(L, 1), &nframes));
  lua_pushinteger(L, (lua_Integer)nframes);
  return 1;
}

static int sys_chill(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cutoff = luaL_checknumber(L, 2);
  int type = (int)luaL_checkinteger(L, 3);
  int graph = (int)luaL_checkinteger(L, 4);
  int rule = (int)luaL_checkinteger(L, 5);
  int n = seams_system_count(box->sys);
  int counts[9];
  int *labels;
  seams_status status;
  if (n <= 0) {
    return luaL_error(L, "system has no frame");
  }
  labels = (int *)calloc((size_t)n, sizeof(int));
  if (labels == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_chill(box->sys, cutoff, type, graph, rule, labels, (size_t)n, counts);
  if (status != SEAMS_OK) {
    free(labels);
    return raise_status(L, status);
  }
  push_ints(L, labels, n);
  push_ints(L, counts, 9);
  free(labels);
  return 2;
}

static int sys_steinhardt(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cutoff = luaL_checknumber(L, 2);
  int type = (int)luaL_checkinteger(L, 3);
  int order = (int)luaL_checkinteger(L, 4);
  int average = (int)luaL_checkinteger(L, 5);
  int n = seams_system_count(box->sys);
  double *ql;
  double *bar = NULL;
  seams_status status;
  if (n <= 0) {
    return luaL_error(L, "system has no frame");
  }
  ql = (double *)calloc((size_t)n, sizeof(double));
  if (average) {
    bar = (double *)calloc((size_t)n, sizeof(double));
  }
  if (ql == NULL || (average && bar == NULL)) {
    free(ql);
    free(bar);
    return luaL_error(L, "out of memory");
  }
  status = seams_steinhardt(box->sys, cutoff, type, order, ql, (size_t)n, bar,
                            bar == NULL ? 0 : (size_t)n);
  if (status != SEAMS_OK) {
    free(ql);
    free(bar);
    return raise_status(L, status);
  }
  push_doubles(L, ql, n);
  if (bar == NULL) {
    lua_pushnil(L);
  } else {
    push_doubles(L, bar, n);
  }
  free(ql);
  free(bar);
  return 2;
}

static int rings_vertices(lua_State *L) {
  RingsBox *box = rings_box(L, 1);
  int index = (int)luaL_checkinteger(L, 2);
  size_t n_out = 0;
  int *atoms = NULL;
  seams_status status = seams_rings_vertices(box->rings, index, NULL, 0, &n_out);
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  atoms = (int *)calloc(n_out > 0 ? n_out : 1, sizeof(int));
  if (atoms == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_rings_vertices(box->rings, index, atoms, n_out, &n_out);
  if (status != SEAMS_OK) {
    free(atoms);
    return raise_status(L, status);
  }
  push_ints(L, atoms, (int)n_out);
  free(atoms);
  return 1;
}

static CageBox *cage_box(lua_State *L, int index) {
  CageBox *box = (CageBox *)luaL_checkudata(L, index, "dseams.Cages");
  if (box->cages == NULL) {
    luaL_error(L, "cages are closed");
  }
  return box;
}

static int cages_close(lua_State *L) {
  CageBox *box = (CageBox *)luaL_checkudata(L, 1, "dseams.Cages");
  seams_cages_destroy(box->cages);
  box->cages = NULL;
  return 0;
}

static void push_cages(lua_State *L, seams_cage_set *cages) {
  CageBox *box = (CageBox *)lua_newuserdata(L, sizeof(CageBox));
  box->cages = cages;
  luaL_getmetatable(L, "dseams.Cages");
  lua_setmetatable(L, -2);
}

static int sys_cages(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  seams_cage_set *cages = NULL;
  seams_status status =
      seams_cages(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                  (int)luaL_checkinteger(L, 4), (int)luaL_checkinteger(L, 5),
                  (int)luaL_checkinteger(L, 6), &cages);
  check_ok(L, status);
  push_cages(L, cages);
  return 1;
}

static int cages_counts(lua_State *L) {
  CageBox *box = cage_box(L, 1);
  int n_hex = 0;
  int n_cubic = 0;
  int n_water = 0;
  int n_both = 0;
  check_ok(L, seams_cages_counts(box->cages, &n_hex, &n_cubic, &n_water, &n_both));
  lua_pushinteger(L, n_hex);
  lua_pushinteger(L, n_cubic);
  lua_pushinteger(L, n_water);
  lua_pushinteger(L, n_both);
  return 4;
}

static int cages_labels(lua_State *L) {
  CageBox *box = cage_box(L, 1);
  int n = (int)luaL_checkinteger(L, 2);
  int *labels;
  seams_status status;
  if (n < 0) {
    return luaL_error(L, "label length is negative");
  }
  labels = (int *)calloc((size_t)(n > 0 ? n : 1), sizeof(int));
  if (labels == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_cages_labels(box->cages, labels, (size_t)n);
  if (status != SEAMS_OK) {
    free(labels);
    return raise_status(L, status);
  }
  push_ints(L, labels, n);
  free(labels);
  return 1;
}

static int sys_signature(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  size_t n_cages = 0;
  size_t n_atoms = 0;
  check_ok(L, seams_signature_cages(box->sys, luaL_checknumber(L, 2),
                                    (int)luaL_checkinteger(L, 3),
                                    (int)luaL_checkinteger(L, 4), luaL_checkstring(L, 5),
                                    &n_cages, &n_atoms));
  lua_pushinteger(L, (lua_Integer)n_cages);
  lua_pushinteger(L, (lua_Integer)n_atoms);
  return 2;
}

static int sys_knn(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  seams_neighbours *nb = NULL;
  seams_status status =
      seams_neighbours_knn(box->sys, (int)luaL_checkinteger(L, 2), luaL_checknumber(L, 3),
                           (int)luaL_checkinteger(L, 4), (int)luaL_checkinteger(L, 5), &nb);
  check_ok(L, status);
  push_neighbours(L, nb);
  return 1;
}

static int sys_rdf(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double rmax = luaL_checknumber(L, 2);
  int bins = (int)luaL_checkinteger(L, 3);
  int type_i = (int)luaL_checkinteger(L, 4);
  int type_j = (int)luaL_checkinteger(L, 5);
  size_t n_out = 0;
  double volume = 0.0;
  double *r = NULL;
  double *g = NULL;
  seams_status status =
      seams_rdf(box->sys, rmax, bins, type_i, type_j, NULL, NULL, NULL, 0, &n_out, &volume);
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  r = (double *)calloc(n_out > 0 ? n_out : 1, sizeof(double));
  g = (double *)calloc(n_out > 0 ? n_out : 1, sizeof(double));
  if (r == NULL || g == NULL) {
    free(r);
    free(g);
    return luaL_error(L, "out of memory");
  }
  status = seams_rdf(box->sys, rmax, bins, type_i, type_j, r, g, NULL, n_out, &n_out, NULL);
  if (status != SEAMS_OK) {
    free(r);
    free(g);
    return raise_status(L, status);
  }
  push_doubles(L, r, (int)n_out);
  push_doubles(L, g, (int)n_out);
  lua_pushnumber(L, volume);
  free(r);
  free(g);
  return 3;
}

static int sys_running_cn(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  size_t n_out = 0;
  double *r = NULL;
  double *cn = NULL;
  seams_status status;
  status = seams_running_cn(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                            (int)luaL_checkinteger(L, 4), (int)luaL_checkinteger(L, 5), NULL,
                            NULL, 0, &n_out);
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  r = (double *)calloc(n_out > 0 ? n_out : 1, sizeof(double));
  cn = (double *)calloc(n_out > 0 ? n_out : 1, sizeof(double));
  if (r == NULL || cn == NULL) {
    free(r);
    free(cn);
    return luaL_error(L, "out of memory");
  }
  status = seams_running_cn(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                            (int)luaL_checkinteger(L, 4), (int)luaL_checkinteger(L, 5), r, cn,
                            n_out, &n_out);
  if (status != SEAMS_OK) {
    free(r);
    free(cn);
    return raise_status(L, status);
  }
  push_doubles(L, r, (int)n_out);
  push_doubles(L, cn, (int)n_out);
  free(r);
  free(cn);
  return 2;
}

static int sys_coordination(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  double cn = 0.0;
  check_ok(L, seams_coordination(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                                 (int)luaL_checkinteger(L, 4), (int)luaL_checkinteger(L, 5),
                                 &cn));
  lua_pushnumber(L, cn);
  return 1;
}

static int sys_hbonds(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  int n_bonds = 0;
  check_ok(L, seams_hbonds(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                           (int)luaL_checkinteger(L, 4), luaL_checknumber(L, 5),
                           luaL_checknumber(L, 6), &n_bonds));
  lua_pushinteger(L, n_bonds);
  return 1;
}

static int sys_density(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  int bins = (int)luaL_checkinteger(L, 3);
  size_t n_out = 0;
  double *z = NULL;
  double *rho = NULL;
  seams_status status =
      seams_density_axis(box->sys, (int)luaL_checkinteger(L, 2), bins,
                         (int)luaL_checkinteger(L, 4), NULL, NULL, 0, &n_out);
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  z = (double *)calloc(n_out > 0 ? n_out : 1, sizeof(double));
  rho = (double *)calloc(n_out > 0 ? n_out : 1, sizeof(double));
  if (z == NULL || rho == NULL) {
    free(z);
    free(rho);
    return luaL_error(L, "out of memory");
  }
  status = seams_density_axis(box->sys, (int)luaL_checkinteger(L, 2), bins,
                              (int)luaL_checkinteger(L, 4), z, rho, n_out, &n_out);
  if (status != SEAMS_OK) {
    free(z);
    free(rho);
    return raise_status(L, status);
  }
  push_doubles(L, z, (int)n_out);
  push_doubles(L, rho, (int)n_out);
  free(z);
  free(rho);
  return 2;
}

static int sys_fingerprint(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  seams_topo *fp = NULL;
  size_t key_len = 0;
  char *key = NULL;
  int counts[13];
  seams_status status =
      seams_fingerprint(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                        (int)luaL_checkinteger(L, 4), (int)luaL_checkinteger(L, 5),
                        (int)luaL_checkinteger(L, 6), (int)luaL_checkinteger(L, 7), &fp);
  if (status != SEAMS_OK) {
    return raise_status(L, status);
  }
  status = seams_fingerprint_key(fp, NULL, 0, &key_len);
  if (status != SEAMS_OK) {
    seams_fingerprint_destroy(fp);
    return raise_status(L, status);
  }
  key = (char *)calloc(key_len + 1, 1);
  if (key == NULL) {
    seams_fingerprint_destroy(fp);
    return luaL_error(L, "out of memory");
  }
  status = seams_fingerprint_key(fp, key, key_len + 1, &key_len);
  if (status == SEAMS_OK) {
    status = seams_fingerprint_rings(fp, counts, 13);
  }
  seams_fingerprint_destroy(fp);
  if (status != SEAMS_OK) {
    free(key);
    return raise_status(L, status);
  }
  lua_pushlstring(L, key, key_len);
  push_ints(L, counts, 13);
  free(key);
  return 2;
}

static int sys_pairs(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  size_t n_pairs = 0;
  int n_cation = 0;
  int n_anion = 0;
  check_ok(L, seams_contact_pairs(box->sys, luaL_checkstring(L, 2), &n_pairs, &n_cation,
                                  &n_anion));
  lua_pushinteger(L, (lua_Integer)n_pairs);
  lua_pushinteger(L, n_cation);
  lua_pushinteger(L, n_anion);
  return 3;
}

static int sys_domains(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  int n_subset = 0;
  int largest = 0;
  double percolation = 0.0;
  check_ok(L, seams_domains(box->sys, luaL_checknumber(L, 2), luaL_checkstring(L, 3),
                            (int)luaL_checkinteger(L, 4), &n_subset, &largest, &percolation));
  lua_pushinteger(L, n_subset);
  lua_pushinteger(L, largest);
  lua_pushnumber(L, percolation);
  return 3;
}

static int sys_prisms(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  int perfect[7];
  int imperfect[7];
  check_ok(L, seams_prisms(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                           perfect, imperfect, 7));
  push_ints(L, perfect, 7);
  push_ints(L, imperfect, 7);
  return 2;
}

static int sys_monolayer(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  int depth = (int)luaL_checkinteger(L, 4);
  int *counts;
  seams_status status;
  if (depth < 3) {
    return luaL_error(L, "max_depth must be at least 3");
  }
  counts = (int *)calloc((size_t)depth + 1, sizeof(int));
  if (counts == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_monolayer(box->sys, luaL_checknumber(L, 2), (int)luaL_checkinteger(L, 3),
                           depth, luaL_checknumber(L, 5), counts, (size_t)depth + 1);
  if (status != SEAMS_OK) {
    free(counts);
    return raise_status(L, status);
  }
  push_ints(L, counts, depth + 1);
  free(counts);
  return 1;
}

static int sys_sites(lua_State *L) {
  SystemBox *box = system_box(L, 1);
  int n = seams_system_count(box->sys);
  int *kinds;
  seams_status status;
  if (n <= 0) {
    return luaL_error(L, "system has no frame");
  }
  kinds = (int *)calloc((size_t)n, sizeof(int));
  if (kinds == NULL) {
    return luaL_error(L, "out of memory");
  }
  status = seams_site_table(box->sys, luaL_checkstring(L, 2), kinds, (size_t)n);
  if (status != SEAMS_OK) {
    free(kinds);
    return raise_status(L, status);
  }
  push_ints(L, kinds, n);
  free(kinds);
  return 1;
}

static void register_type(lua_State *L, const char *name, const luaL_Reg *methods) {
  luaL_newmetatable(L, name);
  lua_pushvalue(L, -1);
  lua_setfield(L, -2, "__index");
  luaL_setfuncs(L, methods, 0);
  lua_pop(L, 1);
}

static const luaL_Reg system_methods[] = {
    {"set_frame", sys_set_frame},
    {"set_ortho", sys_set_ortho},
    {"count", sys_count},
    {"chill_plus", sys_chill_plus},
    {"q4q6", sys_q4q6},
    {"neighbours", sys_neighbours},
    {"rings", sys_rings},
    {"read", sys_read},
    {"chill", sys_chill},
    {"steinhardt", sys_steinhardt},
    {"cages", sys_cages},
    {"signature", sys_signature},
    {"knn", sys_knn},
    {"rdf", sys_rdf},
    {"running_cn", sys_running_cn},
    {"coordination", sys_coordination},
    {"hbonds", sys_hbonds},
    {"density", sys_density},
    {"fingerprint", sys_fingerprint},
    {"contact_pairs", sys_pairs},
    {"domains", sys_domains},
    {"prisms", sys_prisms},
    {"monolayer", sys_monolayer},
    {"site_table", sys_sites},
    {"close", sys_close},
    {"__gc", sys_close},
    {"__close", sys_close},
    {NULL, NULL},
};

static const luaL_Reg neighbour_methods[] = {
    {"count", nb_count},
    {"csr", nb_csr},
    {"close", nb_close},
    {"__gc", nb_close},
    {"__close", nb_close},
    {NULL, NULL},
};

static const luaL_Reg rings_methods[] = {
    {"count", rings_count},
    {"histogram", rings_histogram},
    {"labels", rings_labels},
    {"fraction", rings_fraction},
    {"vertices", rings_vertices},
    {"close", rings_close},
    {"__gc", rings_close},
    {"__close", rings_close},
    {NULL, NULL},
};

static const luaL_Reg module_funcs[] = {
    {"system", l_system},
    {"version", l_version},
    {"abi", l_abi},
    {"status_string", l_status_string},
    {"label_name", l_label_name},
    {"ring_symmetry_name", l_ring_name},
    {"set_num_threads", l_set_threads},
    {"num_threads", l_threads},
    {"last_error", l_last_error},
    {"nframes", l_nframes},
    {NULL, NULL},
};

static const luaL_Reg cage_methods[] = {
    {"counts", cages_counts},
    {"labels", cages_labels},
    {"close", cages_close},
    {"__gc", cages_close},
    {"__close", cages_close},
    {NULL, NULL},
};

static void push_const(lua_State *L, const char *name, int value) {
  lua_pushinteger(L, value);
  lua_setfield(L, -2, name);
}

int luaopen_dseams_core(lua_State *L) {
  register_type(L, "dseams.System", system_methods);
  register_type(L, "dseams.Neighbours", neighbour_methods);
  register_type(L, "dseams.Rings", rings_methods);
  register_type(L, "dseams.Cages", cage_methods);
  luaL_newlib(L, module_funcs);
  push_const(L, "PUBLISHED", SEAMS_CHILL_PUBLISHED);
  push_const(L, "OVITO", SEAMS_CHILL_OVITO);
  push_const(L, "CUTOFF", SEAMS_GRAPH_CUTOFF);
  push_const(L, "MUTUAL", SEAMS_GRAPH_MUTUAL);
  push_const(L, "SEEDED", SEAMS_GRAPH_SEEDED);
  push_const(L, "WATER_O", SEAMS_SITE_WATER_O);
  push_const(L, "OK", SEAMS_OK);
  push_const(L, "ABI", SEAMS_ABI);
  return 1;
}
