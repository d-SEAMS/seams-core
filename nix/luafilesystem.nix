# [[file:~/Git/Github/C++/Forks/structureFactor/literateNix.org::luafilesystem][luafilesystem]]
fileSystem = buildLuaPackage {
  name = "filesystem-1.6.2";
  src = fetchurl {
    url = "https://github.com/keplerproject/luafilesystem/archive/v1_6_2.tar.gz";
    sha256 = "1n8qdwa20ypbrny99vhkmx8q04zd2jjycdb5196xdhgvqzk10abz";
  };
  meta = {
    homepage = "https://github.com/keplerproject/luafilesystem";
    hydraPlatforms = stdenv.lib.platforms.linux;
    maintainers = with maintainers; [ flosse ];
  };
};
# luafilesystem ends here
