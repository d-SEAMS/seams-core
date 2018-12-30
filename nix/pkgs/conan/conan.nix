# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::*Conan][Conan:3]]
{ lib, python3, fetchpatch, git }:

let newPython = python3.override {
  packageOverrides = self: super: {
    distro = super.distro.overridePythonAttrs (oldAttrs: rec {
      version = "1.2.0";
      src = oldAttrs.src.override {
        inherit version;
        sha256 = "1vn1db2akw98ybnpns92qi11v94hydwp130s8753k6ikby95883j";
      };
    });
    node-semver = super.node-semver.overridePythonAttrs (oldAttrs: rec {
      version = "0.2.0";
      src = oldAttrs.src.override {
        inherit version;
        sha256 = "1080pdxrvnkr8i7b7bk0dfx6cwrkkzzfaranl7207q6rdybzqay3";
      };
    });
    future = super.future.overridePythonAttrs (oldAttrs: rec {
      version = "0.16.0";
      src = oldAttrs.src.override {
        inherit version;
        sha256 = "1nzy1k4m9966sikp0qka7lirh8sqrsyainyf8rk97db7nwdfv773";
      };
    });
    tqdm = super.tqdm.overridePythonAttrs (oldAttrs: rec {
      version = "4.28.1";
      src = oldAttrs.src.override {
        inherit version;
        sha256 = "1fyybgbmlr8ms32j7h76hz5g9xc6nf0644mwhc40a0s5k14makav";
      };
    });
  };
};

in newPython.pkgs.buildPythonApplication rec {
  version = "1.9.1";
  pname = "conan";

  src = newPython.pkgs.fetchPypi {
    inherit pname version;
    sha256 = "0mn69ps84w8kq76zba2gnlqlp855a6ksbl1l6pd1gkjlp9ry0hnf";
  };
  checkInputs = [
    git
  ] ++ (with newPython.pkgs; [
    nose
    parameterized
    mock
    webtest
    codecov
  ]);

  propagatedBuildInputs = with newPython.pkgs; [
    requests fasteners pyyaml pyjwt colorama patch
    bottle pluginbase six distro pylint node-semver
    future pygments mccabe deprecation tqdm
  ];

  checkPhase = ''
    export HOME="$TMP/conan-home"
    mkdir -p "$HOME"
  '';

  meta = with lib; {
    homepage = https://conan.io;
    description = "Decentralized and portable C/C++ package manager";
    license = licenses.mit;
    platforms = platforms.linux;
  };
}
# Conan:3 ends here
