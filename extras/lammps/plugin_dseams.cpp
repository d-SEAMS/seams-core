#include "compute_dseams.h"
#include "lammps.h"
#include "lammpsplugin.h"
#include "version.h"

using namespace LAMMPS_NS;

static void *compute_dseams_creator(LAMMPS *lmp, int argc, char **argv) {
  return static_cast<void *>(new ComputeDseams(lmp, argc, argv));
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc) {
  lammpsplugin_t plugin{};
  auto register_plugin = reinterpret_cast<lammpsplugin_regfunc>(regfunc);
  plugin.version = LAMMPS_VERSION;
  plugin.style = "compute";
  plugin.name = "dseams";
  plugin.info = "d-SEAMS CHILL+ per-atom compute";
  plugin.author = "d-SEAMS";
  plugin.creator.v2 =
      reinterpret_cast<lammpsplugin_factory2 *>(&compute_dseams_creator);
  plugin.handle = handle;
  register_plugin(&plugin, lmp);
}
