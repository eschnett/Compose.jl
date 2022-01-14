# Compose.jl

Julia utilities for [CompOSE](https://compose.obspm.fr/home) CompStar
Online Supernovae Equations of State. CompOSE provides tabulated
equations of state. This package helps prepare these tables for
hydrodynamics simulations. Currently, they find and correct
inconsisitencies in these tables.

# `dE/dT`

The tables are indexed by baryon number density `n_b`, temperature
`T`, and hadronic charge fraction `Y_q`. One of the table entries is
the internal energy per baryon `E`. To be thermodynamically
consistent, the condition `dE/dT > 0` needs to hold.

This condition is also important when converting between conserved and
primitive variables in a hydrodynamics simulation. If it doesn't hold,
then this conversion is not well defined any more.

Unfortunately, the tables might contain outliers where this condition
doesn't hold, especially for extreme values of `n_b`. This package
finds these outliers and replaces them by interpolated values,
ensuring that this condition holds everywhere.

Here is some example output:
```
[ Info: Read EOS table with (164, 163, 51) entries for 8 quantities
[ Info: Looking for outliers...
[ Info: dE/dT ≤ 0 at 306 of 1363332 points (0.022445009726170882%)
[ Info: Correcting outliers...
[ Info: Corrected all outliers
[ Info: Writing result...
[ Info: Done.
```
