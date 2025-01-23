```@meta
CurrentModule = CherenkovMediumBase
```

# CherenkovMediumBase

Welcome to the documentation for [CherenkovMediumBase](https://github.com/PLEnuM-group/CherenkovMediumBase.jl).

This package provides a set of tools and abstractions for working with medium properties in physics simulations.

## Installation

To install the package, use the following command:

```julia
using Pkg
Pkg.add("CherenkovMediumBase")
```

## Usage

Here is a simple example of how to use the package:

```julia
using CherenkovMediumBase

# Define a custom medium
struct CustomMedium <: MediumProperties end

# Implement required methods
function CherenkovMediumBase.pressure(medium::CustomMedium)
    return 1.0 # Example pressure
end

# ...implement other required methods...

# Create an instance of the custom medium
medium = CustomMedium()

# Use the medium in calculations
println(CherenkovMediumBase.pressure(medium))
```

## Documentation

The following sections provide detailed documentation for the package:

```@index
```