{ clangStdenv, fetchurl }:

clangStdenv.mkDerivation rec {
  name = "rang-${version}";
  version = "3.1.0";

  src = fetchurl {
    url = "https://github.com/agauniyal/rang/releases/download/v${version}/rang.hpp";
      sha256 = "190hvbpcvnzwiaiahcn6cypc0n1q65cy974q10yx9zx8343i11z4";
  };

  # This is a header only library. No unpacking needed. Seems like we need to create
  # _some_ folder, otherwise we get errors.
  unpackCmd = "mkdir fake_dir";

  installPhase = ''
    mkdir -p $out/include/
    cp ${src} $out/include/rang.hpp
  '';

  meta = {
    description = "
A Minimal, Header only Modern c++ library for terminals";
    homepage = https://agauniyal.github.io/rang/;
  };
}
