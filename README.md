# AbstractMediumProperties

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chrhck.github.io/AbstractMediumProperties.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://chrhck.github.io/AbstractMediumProperties.jl/dev/)
[![Build Status](https://github.com/chrhck/AbstractMediumProperties.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/chrhck/AbstractMediumProperties.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/chrhck/AbstractMediumProperties.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/chrhck/AbstractMediumProperties.jl)

AbstractMediumProperties provides abstract types and functions for calculating various properties of different media. It provides functionalities to compute properties such as refractive index, scattering length, absorption length, and more, based on the medium's characteristics.

## Main Features

- Calculation of phase and group refractive indices
- Scattering and absorption length computations
- Cherenkov angle and group velocity calculations
- Support for different scattering functions including Henyey-Greenstein and Simplified-Liu models

## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add("AbstractMediumProperties")
```

## Implementing a Concrete MediumProperties Subtype

To implement a concrete subtype of `MediumProperties`, you need to define a new type and implement the required methods. Here is an example:

```julia
using AbstractMediumProperties

struct Water <: MediumProperties
    salinity::Float64
    temperature::Float64
    pressure::Float64
end

function phase_refractive_index(wavelength, medium::Water)
    return refractive_index_fry(wavelength; salinity=medium.salinity, temperature=medium.temperature, pressure=medium.pressure)
end

function dispersion(wavelength, medium::Water)
    return dispersion_fry(wavelength; salinity=medium.salinity, temperature=medium.temperature, pressure=medium.pressure)
end

# Implement other required methods similarly
```

This example defines a `Water` type and implements the `phase_refractive_index` and `dispersion` methods for it. You can implement other required methods similarly.

## Usage

Here is a simple example of how to use the package:

```julia
using AbstractMediumProperties

# Define a medium (example)
medium = Water()

# Calculate phase refractive index at a given wavelength
wavelength = 500.0 # nm
n_phase = phase_refractive_index(wavelength, medium)

println("Phase refractive index at $wavelength nm: $n_phase")
```

For more detailed usage and examples, please refer to the [documentation](https://chrhck.github.io/AbstractMediumProperties.jl/stable/).



