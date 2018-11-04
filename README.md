# ARPACK for Visual Studio

This is a Visual Studio solution for [ARPACK](https://github.com/opencollab/arpack-ng).

## Instructions

The repository contains a F2C'ed version of ARPACK. It was obtained from https://github.com/opencollab/arpack-ng. Additionally, ARPACK++ headers are needed to compile the project. Download the latest version from https://github.com/m-reuter/arpackpp and extract it to the `src/ARPACK` folder (a subfolder `arpackpp` should be created). To check if everything is in its right place, make sure that the file `src/ARPACK/arpackpp/README.md` exists.

Replace `ARPACK\arpackpp\include\arerror.h` with `ARPACK\arerror.h`.

Pre-compiled binaries for windows users can be found [here](http://wo80.bplaced.net/math/packages.html).

## Thread Safety

At the moment, the ARPACK C code created by F2C is **NOT** thread-safe!

## Why?

The project was created to maintain ARPACK builds matching the [CSparse.Interop](https://github.com/wo80/csparse-interop) bindings for C#.
