# ARPACK for Visual Studio

This is a Visual Studio solution for [ARPACK](https://github.com/opencollab/arpack-ng).

## Instructions

The repository contains a f2c'ed version of ARPACK. It was be obtained from https://github.com/opencollab/arpack-ng. Additionally, ARPACK++ headers are needed to compile the project. Download the latest version from https://github.com/m-reuter/arpackpp and extract it to the `src/ARPACK` folder (a subfolder `arpackpp` should be created). To check if everything is in its right place, make sure that the file `src/ARPACK/arpackpp/README.md` exists.

## Why?

The project was created to maintain ARPACK builds matching the [CSparse.Interop](https://github.com/wo80/csparse-interop) bindings for C#.
