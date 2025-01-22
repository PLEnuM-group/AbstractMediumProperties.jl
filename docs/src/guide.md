```@meta
CurrentModule = AbstractMediumProperties
```

# User Guide

This guide provides an overview of how to use the `AbstractMediumProperties` package.

## Defining a Custom Medium

To define a custom medium, create a new type that inherits from `MediumProperties` and implement the required methods.

```julia
using AbstractMediumProperties

struct CustomMedium <: MediumProperties end

function AbstractMediumProperties.pressure(medium::CustomMedium)
    return 1.0 # Example pressure
end

# Implement other required methods...

medium = CustomMedium()
println(AbstractMediumProperties.pressure(medium))
```

## Advanced Usage

For advanced usage, refer to the API documentation and explore the available functions and types.

```@index
```