Introduction
============
`SOFT <http://ft.nephy.chalmers.se/~hoppe/soft/>`_ is a synthetic radiation
diagnostic designed to be used to study synchrotron radiation and bremsstrahlung
from runaway electrons in tokamaks. This documentation is for the rewrite of
SOFT, called `SOFT2` (or `SOFTv2`).

Improvements in SOFTv2 from the older version include:

- Up to 100x faster simulation times
- Improved models of the radiation
- More robust implementations of key algorithms
- Wide selection of analytical magnetic fields (though support for numerical magnetic fields of course remains)
- CODE/NORSE momentum-space distributions can be taken as input directly
- Ability to "include" configuration files in other configuration files
- Simultaneous support for HDF5 and MAT files

The many new features however also comes at a cost: the input format of
the old version of SOFT is no longer supported. All old config files thus have
to be converted (manually) to the new format.

SOFTv2 is written entirely in C++.
