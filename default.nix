# default.nix
with import <nixpkgs> {};
stdenv.mkDerivation {
    name = "mpi_rust"; # Probably put a more meaningful name here
    buildInputs = [clang
    llvmPackages.libclang.lib
    cfitsio
    pkg-config
    necpp

    autoconf
    automake
    libtool
    cmake
    xorg.libX11
    xorg.libXrandr
    xorg.libXinerama
    xorg.libXcursor
    xorg.libXxf86vm
    xorg.libXi
    libGL
    libGL.out
    libGLU
    libGLU.out
    freeglut
    freeglut.out
    libsForQt5.qt5.qtwayland
    ];
    hardeningDisable = [ "all" ];
    #buildInputs = [gcc-unwrapped gcc-unwrapped.out gcc-unwrapped.lib];
    LIBCLANG_PATH = llvmPackages.libclang.lib+"/lib";
    LD_LIBRARY_PATH= libGL+"/lib";
    QT_QPA_PLATFORM="wayland";
}
