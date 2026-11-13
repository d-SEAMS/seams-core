/* libFuzzer harness for seams::cfg dotenv parsing. */
#include <seams_config.hpp>

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t *data, std::size_t size) {
  const std::string path = "fuzz-seams-config.env";
  {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
      return 0;
    }
    out.write(reinterpret_cast<const char *>(data), static_cast<std::streamsize>(size));
  }
  try {
    seams::cfg::applyFile(path, true);
    (void)seams::cfg::load(path);
  } catch (const std::invalid_argument &) {
  } catch (const std::runtime_error &) {
  }
  std::remove(path.c_str());
  return 0;
}
