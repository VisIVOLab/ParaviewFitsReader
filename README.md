# FitsReader Plugin
ParaView/VTK Reader for visualization of FITS format.

The plugin requires a ParaView build compiled with MPI support and [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) installed.

The CMake variable `CFITSIO_ROOT_DIR` can be used to select an installation path for CFITSIO if it is not detected automatically.
