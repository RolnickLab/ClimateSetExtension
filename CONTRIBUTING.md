# Contributing to this repository

## Design patterns
Two main considerations should be made when contributing to this package.

First, a polymorphic approach, using abstract classes and their concrete implementation,
should be prioritized in order to increase maintainability and extensibility.

Therefore, new additions should try to follow this design pattern and either implement
new concrete classes or create new abstract classes and their implementations for 
completely new behavior or needs.

Secondly, a dependency-injection approach is to be preferred, as well as a composition 
approach when creating new modules or extending existing ones.

## Tests

New contributions should include appropriate tests.

## Docstring and type hinting

Docstring format should follow the Numpy standard and type hinting should be used
as per the PEP8 standard : https://docs.python.org/3/library/typing.html