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
  const int *offsets = NULL;
  const int *cols = NULL;
  int n_offsets = 0;
  int nnz = 0;
  seams_status status = seams_neighbours_csr(box->nb, &offsets, &n_offsets, &cols, &nnz);
  check_ok(L, status);
  push_ints(L, offsets, n_offsets);
  if (cols == NULL || nnz == 0) {
    lua_createtable(L, 0, 0);
  } else {
    push_ints(L, cols, nnz);
  }
  lua_pushinteger(L, nnz);
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
  luaL_newlib(L, module_funcs);
  push_const(L, "PUBLISHED", SEAMS_CHILL_PUBLISHED);
  push_const(L, "OVITO", SEAMS_CHILL_OVITO);
  push_const(L, "CUTOFF", SEAMS_GRAPH_CUTOFF);
  push_const(L, "MUTUAL", SEAMS_GRAPH_MUTUAL);
  push_const(L, "OK", SEAMS_OK);
  push_const(L, "ABI", SEAMS_ABI);
  return 1;
}
