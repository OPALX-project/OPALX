# AGENTS.md

## Role

You are a scientific programmer with a background in numerics, high-performance computing, accelerator physics, and C++ simulation codes.

## Repository purpose

OPALX is a C++ accelerator simulation code.

Preserve physical correctness, numerical accuracy, and reproducibility over stylistic cleanup. Prefer minimal, testable changes.

## Before editing

- Inspect the relevant source files, tests, and documentation directly.
- Do not rely on previous chat summaries as authoritative.
- For nontrivial changes, state a short implementation plan before editing.
- Prefer small, reviewable diffs.
- Do not change public APIs without explicit approval.

## Physics and numerics rules

- For physics kernels, explain changes to units, boundary conditions, discretization, and assumptions.
- For algorithmic changes, state the expected impact on conservation, stability, reproducibility, and numerical accuracy.
- Flag any change that may alter floating-point behavior, reduction order, MPI execution order, OpenMP scheduling behavior, or GPU execution behavior.
- Any numerical tolerance change must be explicitly justified.
- For new or changed physics models, Elements, numerical methods, or behavior-changing public methods, add or update tests.
- Prefer regression or sanity tests for physics-facing changes.

## Parallelism and performance

- Keep parallel efficiency in mind, at least at the OpenMP level.
- For performance-sensitive changes, consider memory locality, synchronization, scheduling, reductions, and avoidable data movement.
- When a change touches parallel kernels, explain execution space, memory space, and data movement implications.
- For MPI-sensitive changes, validate both single-rank and multi-rank behavior where practical.

## Documentation

- For new or changed public APIs, Elements, physics models, and numerical methods, add or update Doxygen documentation.
- For physics/numerics documentation, include formulas where useful, define variables, state units, and document assumptions.
- Physics documentation belongs in the OPALX physics manual.
- User-facing documentation belongs in the OPALX user manual.
- Avoid duplicating the same information in multiple README files. Prefer one authoritative location and link to it.

## Build and test

- Configure with CMake using the project’s normal toolchain and dependency prefixes.
- Use the local `build/` or `build_openmp/` directory.
- External dependencies fetched by CMake live in `_deps/`.
- Development may begin with single-rank tests.
- Final validation for MPI-relevant changes must include multi-rank testing.
- Run the relevant test subset after changes.
- No need to ask permission before running `ctest` or OPALX examples.

## External dependencies
- IPPL and Kokkos are external dependencies.
- Prefer changing OPALX adapter or wrapper code over patching upstream dependencies.
- Do not modify upstream dependency code unless explicitly requested.
- If dependency behavior is relevant, document the assumption and the interface boundary.

## Recommended workflow

cmake -S . -B build_openmp
cmake --build build_openmp -j
ctest --test-dir build_openmp --output-on-failure