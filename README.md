# COESA

| **Build Status** | **Code Coverage** |
| :---: | :---: |
| [![Build Status](https://travis-ci.org/danielmatz/COESA.jl.svg?branch=master)](https://travis-ci.org/danielmatz/COESA.jl) [![Build status](https://ci.appveyor.com/api/projects/status/lw6emr2km0m7whgx/branch/master?svg=true)](https://ci.appveyor.com/project/danielmatz/coesa-jl/branch/master) | [![codecov.io](http://codecov.io/github/danielmatz/COESA.jl/coverage.svg?branch=master)](http://codecov.io/github/danielmatz/COESA.jl?branch=master) |

The U.S. Committee on Extension to the Standard Atmosphere (COESA) atmosphere
model, also known as the U.S. Standard Atmosphere, 1976.

## Usage

The main function is `atmosphere(z)`, which returns the atmospheric state at
geometric altitude `z`, where `z` has units of meters.

The atmospheric state object that is returned has the following accessors:
`density`, `pressure`, `temperature`, `speed_of_sound`, and
`mean_molecular_weight`.

```julia
using COESA
atmos = atmosphere(123.0)
rho = density(atmos)
T = temperature(atmos)
P = pressure(atmos)
c = speed_of_sound(atmos)
M = mean_molecular_weight(atmos)
```

## Implementation

For altitudes below 86km, the equations from the original report are used.  The
published tables have a discontinuity in the mean molecular weight and
temperature at 86km.  The report outlines how to blend out the discontinuity.
We implement this blending here.  In this region, the model output matches the
published tables nearly perfectly.  Only a few of the tested altitudes show
differences.  It is only in the pressure values, and in these cases we only miss
the least significant digit.

For the region above 86km, the temperature is computed using the equations from
the original report, but the equations for the mean molecular weight and
pressure are much more complex.  Instead, we adapt the method from Regan's
_Re-Entry Vehicle Dynamics_ to interpolate the tabulated mean molecular weight
and pressure.  We use a quadratic interpolation on the mean molecular weight and
the natural log of the pressure.  Because we are using interpolation, the
computed pressure and density occasionally don't match the published values.
The values are very small, and we still have at least 2 significant digits.

The density is computed as outlined by the original report.

The speed of sound is computed as outlined by the original report for altitudes
below 86km.  For altitudes above 86km, the speed of sound at 86km is used.

## References

1. _U.S. Standard Atmosphere, 1976_. Stock No. 003-017-00323-0.
http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539.pdf.
2. Regan, F.J., _Re-Entry Vehicle Dynamics_, AIAA Press, New York, 1984.
