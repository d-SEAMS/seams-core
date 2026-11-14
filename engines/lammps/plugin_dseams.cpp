// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// `plugin load dseamsplugin.so` registers compute dseams in a stock LAMMPS
// built with the PLUGIN package.

#include "compute_dseams.h"
#include "lammps.h"

#if __has_include("lammpsplugin.h")
#include "lammpsplugin.h"
#else
// Installed LAMMPS trees do not always ship lammpsplugin.h. This is the
// plugin record that PLUGIN reads, unchanged since the 2021 releases.
extern "C" {
typedef void *(lammpsplugin_factory1)(void *);
typedef void *(lammpsplugin_factory2)(void *, int, char **);
typedef struct {
  const char *version;
  const char *style;
  const char *name;
  const char *info;
  const char *author;
  union {
    lammpsplugin_factory1 *v1;
    lammpsplugin_factory2 *v2;
  } creator;
  void *handle;
} lammpsplugin_t;
typedef void (*lammpsplugin_regfunc)(lammpsplugin_t *, void *);
}
#endif

using namespace LAMMPS_NS;

static void *compute_dseams_creator(LAMMPS *lmp, int argc, char **argv) {
  return new ComputeDseams(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc) {
  lammpsplugin_t plugin{};
  auto register_plugin = reinterpret_cast<lammpsplugin_regfunc>(regfunc);
  // The loader compares this string with the running LAMMPS.
  plugin.version = static_cast<LAMMPS *>(lmp)->version;
  plugin.style = "compute";
  plugin.name = "dseams";
  plugin.info = "d-SEAMS CHILL+ per-atom labels";
  plugin.author = "d-SEAMS developers";
  plugin.creator.v2 = reinterpret_cast<lammpsplugin_factory2 *>(&compute_dseams_creator);
  plugin.handle = handle;
  register_plugin(&plugin, lmp);
}
