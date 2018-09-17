This directory should contain external dependencies in the corresponding `x86` or `x64` subfolders.

The Visual Studio solution expects the following static libraries to be present:
 * libf2c.lib
 * blas.lib
 * lapack.lib
 * libsuperlu.lib

The SuperLU dependency can be build using [vs-superlu](https://github.com/wo80/vs-superlu).
