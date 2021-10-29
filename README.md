# yPlusStats

**TL;DR**: Do not write a volScalarField every time the yPlus function object is called.

## Requirements

* git
* OpenFOAM v1906+

## Rationale

When calculating y+, OpenFOAM generates two outputs: one containing the min,
max, and average values on each patch for each time the yPlus function is
called; and a second one storing the y+ for the _entire domain_ as a
volScalarField. This code disables the second output, which may be of the
interest for some users.

The problem of writing a y+ volScalarField arises when one uses `purgeWrite` or
a `writeInterval` different from the one in `system/controlDict`. In these cases,
OpenFOAM creates various time folders containing only the y+ volScalarField, which
may lead to the following problems:

* The total number of folders in one's simulation can rapidly increase. In some
  cases, using tiny time steps and a large amount of subdomains can
  exponentially grow the number of folders created. And thus, the simulation
  may reach some sort of limit that prohibits writing to disk.

* If the write interval for the y+ is different from the one for other flow
  variables, restarting the simulation may require the user to manually delete
  some time folders in which only the y+ volScalarField is stored.

## Install

  1. Load the OpenFOAM environment.

  2. Create the user source folder.

    mkdir -p ${FOAM_RUN}/../src

  3. Clone the repository

    cd ${FOAM_RUN}/../src
    git clone https://github.com/gabrielbdsantos/yPlusStats
    cd yPlusStats

  4. Compile the library.

    wmake all

## Disclaimer and License

The code is **not** approved or endorsed by OpenCFD Ltd., the producer of the
OpenFOAM software and owner of the OprenFOAM and OpenCFD trade marks. Any
implementations presented here are bound to the same license as OpenFOAM: the
GNU Public License v3. See the [LICENSE](./LICENSE) for more information.
