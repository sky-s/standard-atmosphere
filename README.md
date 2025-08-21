# Standard Atmosphere Functions
Standard atmosphere gas properties. Support for n-dim inputs, non-standard atmospheres, units, etc.

[![View Standard Atmosphere Functions on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/28135)

Standard atmosphere functions based on the 1976 Standard Atmosphere. Returns density, speed of sound, temperature, pressure, and viscosity for a given altitude input up to 86 km.

Functions are designed to be useful for those designing and analyzing aircraft and have the following features:
- Inputs may be scalar, vector, matrix, or n-dimensional arrays. Functions are vectorized and fast for computing conditions at a large number of points simultaneously. Especially for n-dimensional problems, it is faster than the built-in atmosisa that comes with the aerospace toolbox.
- Temperature offset option for non-standard atmospheres, e.g analyzing hot day aircraft performance.
- Absolute temperature option for non-standard atmospheres, e.g for when outside air temperature is known.
- Non-standard atmosphere profiles (hot day, cold day, tropical, polar) via `atmos` function based on SAE AS210.
- Either SI or imperial units (and easy to set your preferred default).
- Units consistency can be enforced by using the Physical Units Toolbox, reducing errors and making code clearer.
- Returns everything needed to easily determine important parameters such as dynamic pressure, Mach number, Reynolds number, stagnation temperature, etc.
- Option for geometric instead of geopotential altitude input.
- Density altitude function allows reverse lookup of altitude based on gas properties.
- Stripped-down, troposphere-only function included for when computation speed is a priority.

## Usage Examples

### Standard vs Non-Standard Atmospheres
```matlab
% Compare atmospheres at 10 km altitude
h = 10000; % meters
[rho_std, a_std, T_std] = atmos(h, 'atmosphere', 'standard');
[rho_hot, a_hot, T_hot] = atmos(h, 'atmosphere', 'hot');
[rho_cold, a_cold, T_cold] = atmos(h, 'atmosphere', 'cold');

fprintf('At 10 km:\n');
fprintf('Standard: T=%.1f K, rho=%.3f kg/m³\n', T_std, rho_std);
fprintf('Hot:      T=%.1f K, rho=%.3f kg/m³\n', T_hot, rho_hot);
fprintf('Cold:     T=%.1f K, rho=%.3f kg/m³\n', T_cold, rho_cold);
```

References: 
- U.S. Standard Atmosphere, 1976
- SAE AS210 - Environmental Conditions for Aerospace  
- ESDU 77022; www.pdas.com/atmos.html
