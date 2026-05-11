# **CompOSE**

The online service **CompOSE** (http://compose.obspm.fr) provides information and
data tables for different equations of states (EoS) ready for further use
in astrophysical applications, nuclear physics and beyond. This service has
three major purposes:

- **CompOSE** is a repository of EoS tables in a common format for direct usage
  with information on a large number of thermodynamic properties, on the
  chemical composition of dense matter and, if available, on microphysical
  quantities of the constituents.
- **CompOSE** allows to interpolate the
  provided tables using different schemes to obtain the relevant quantities,
  selected by the user, for grids that are tailored to specific applications.
- **CompOSE** can provide additional thermodynamic quantities, which
  are not stored in the original data tables, and on further quantities, which
  characterize an EoS such as nuclear matter parameters and compact star
  properties.

The format of the files as well as the calculational mesh is mainly determined
according to the needs of scientific groups performing extensive numerical
simulations of astrophysical objects.

We cannot offer an online service for all features of **CompOSE**, e.g. to run
all codes online.  This is mostly due to limitations in storage and
computation times but also gives better control on avoiding unphysical input
parameters. However, we offer several computational tools that allow the user
to extract the data from the tables that are relevant for her/him. The
present release of the **CompOSE** tools contains the following files

- Makefile
  A Makefile using the gnu make utility which allows to build the program
  compose from the existing source files (see below). At present it employs
  the gfortran and gcc compilers. A flag can be set in the Makefile to
  enable/disable the compilation of the HDF5 output part.

* compose.f90
  Main compose code.

* composemodules.f90
  Contains all the modules used by compose apart from those specific to the
  HDF5 output part.

* hdf5compose.f90
  Interface for HDF5 output.

* hdf5writecompose.c
  HDF5 output routines.

* hdf5readcompose.c
  Example routines for reading the HDF5 output.

**CompOSE** is a service initiated by the CompStar European network and maintained
mainly by C. Ishizuka, T. Klaehn, M. Marcini, M. Oertel, and S. Typel (contact develop@compose.obspm.fr).


### INSTALL :

I the standard case type "make" or "make compose"

## Installation with HDF5
if you want install Compose without HDF5 !

In a terminal :
```bash
export USE_HDF5 = 0
make
```
