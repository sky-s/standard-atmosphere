# Standard Atmosphere Functions
Standard atmosphere gas properties. Support for n-dim inputs, non-standard atmospheres, units, etc.

[![View Standard Atmosphere Functions on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/28135)

Standard atmosphere functions based on the 1976 Standard Atmosphere. Returns density, speed of sound, temperature, pressure, and viscosity for a given altitude input up to 86 km.

Functions are designed to be useful for those designing and analyzing aircraft and have the following features:
- Inputs may be scalar, vector, matrix, or n-dimensional arrays. Functions are vectorized and fast for computing conditions at a large number of points simultaneously. Especially for n-dimensional problems, it is faster than the built-in atmosisa that comes with the aerospace toolbox.
- Temperature offset option for non-standard atmospheres, e.g analyzing hot day aircraft performance.
- Absolute temperature option for non-standard atmospheres, e.g for when outside air temperature is known.
- Either SI or imperial units (and easy to set your preferred default).
- Units consistency can be enforced by using the Physical Units Toolbox, reducing errors and making code clearer.
- Returns everything needed to easily determine important parameters such as dynamic pressure, Mach number, Reynolds number, stagnation temperature, etc.
- Option for geometric instead of geopotential altitude input.
- Density altitude function allows reverse lookup of altitude based on gas properties.
- Stripped-down, troposphere-only function included for when computation speed is a priority.

References: ESDU 77022; www.pdas.com/atmos.html
