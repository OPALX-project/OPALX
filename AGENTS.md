# AGENTS.md

## Repository purpose
OPALX is a C++ accelerator simulation code. Preserve physical correctness over stylistic cleanup.

## Coding
- Always make a Doxygen documentation with possible formulas 
- For each Element/Methods/new Physic add a unit test
- Keep am eye on parallel efficincy at least on the openmp level

## Build
- Configure with CMake using the project’s normal toolchain and dependency prefixes.
- Build in the local build or build_openmp directory.
- Test must be done single and multirank but development can be done first single rank
- build instructions can be found online.
- Run the relevant test subset after changes.
- External dependencies are in the _deps directory fetched with cmake

## External dependencies
- IPPL and Kokkos are external dependencies.
- Prefer changing OPALX adapter/wrapper code over patching upstream dependencies.
- When a change touches parallel kernels, explain execution space, memory space, and data movement implications.

## Physics / numerics rules
- For algorithmic changes, state expected impact on conservation, stability, and reproducibility.
- Flag any change that may alter floating-point behavior or MPI/GPU execution order.
- Add or update at least one regression/sanity test for physics-facing changes.
- Keep am eye on parallel efficincy at least on the openmp level
- All physics need to be documented at git/OPALX-project.github.io/gamma-gamma (follow the expamles)


## Permissions
- you have permission in ~/git and /tmp
- no need to ask for edits, awk and python tasks.
- repository pushes need permission

## Python
- use the git/opalx/.venv-h6
- organize data in pandas dataframes
- use always high quality figures
