# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What FEMFFUSION Is

FEMFFUSION is a C++ finite element method (FEM) code for nuclear reactor neutronics. It solves the multigroup neutron transport equation using diffusion or SPN (Simplified Spherical Harmonics) approximations in 1D, 2D, and 3D. It is built on top of deal.II, PETSc, and SLEPc.

## Build and Run

**Dependencies (Ubuntu 24.04):**
```bash
sudo apt-get install libdeal.ii-dev petsc-dev slepc-dev cmake gmsh
```

**Build:**
```bash
cmake .
make release
make -j$(nproc)
```

**Run a problem:**
```bash
./femffusion.exe -f examples/2D_BIBLIS/biblis_SP3.prm
```

**Run all tests** (must build in Release mode first):
```bash
./femffusion.exe -t
```

**Run all examples:**
```bash
./run_examples.sh
```

**Flags:**
- `-f <file.prm>` — specify input parameter file
- `-t` / `--test` — run the test suite
- `-v` / `--verbose` — verbose output
- `-s` / `--silent` — suppress output
- `-p` / `--performance` — run performance benchmarks

**Parallel core control** (if code uses too many cores unintentionally):
```bash
export OMP_NUM_THREADS=1
# or
mpirun -n 1 ./femffusion.exe -f file.prm
```

## Code Formatting

Uses `.clang-format` (LLVM-based with tabs, indent width 4):
```bash
clang-format -i src/**/*.cc include/**/*.h
```

## Architecture

### Template Pattern

All solver classes are templated on two integers: `dim` (spatial dimension: 1, 2, or 3) and `n_fe_degree` (FEM polynomial degree: 1–5). This means all 15 combinations are compiled explicitly. The `main()` function in `src/femffusion.cc` is a large dispatch tree that reads `Dimension` and `FE_Degree` from the `.prm` file and constructs the right template instantiation.

### Solver Hierarchy

Each problem type has a **static solver** that is constructed first. The static solver computes the steady-state eigenvalue problem (k-eigenvalue or α-eigenvalue). Downstream solvers take a reference to a completed static solver:

```
StaticDiffusion<dim, deg>     — diffusion approximation
StaticSPN<dim, deg>           — SPN approximation
StaticSDPN<dim, deg>          — Simplified Double-PN approximation
StaticFullSPN<dim, deg>       — Full SPN approximation

TimeNeutronDiffusion<dim, deg>(prm, static_prob)  — time-dependent diffusion
TimeNeutronSPN<dim, deg>(prm, static_prob)        — time-dependent SPN
NoiseDiffusion<dim, deg>(prm, static_prob)        — frequency-domain noise
NoiseSPN<dim, deg>(prm, static_prob)              — noise with SPN
NoiseFullSPN<dim, deg>(prm, static_prob)          — noise with Full SPN
ROMKinetics<dim, deg>(prm, static_prob)           — reduced-order model kinetics
ROMStatic<dim, deg>(prm, static_prob)             — reduced-order model static
```

### Key Subdirectories

- `src/eps_solvers/` — eigenvalue problem solvers: `eps_bifpam`, `eps_newton`, `eps_power_it`, `eps_gd`, `eps_solver_2g`, plus preconditioners
- `src/matrix_operators/` — PETSc matrix assembly and matrix-free implementations; separate variants for diffusion, SPN, SDPN, noise (complex), time-dependent, and ROM (LUPOD)
- `src/io/` — cross-section input (`materials.cc`, formats: XS2G, XSEC, XML, Valkin), geometry (`input_geom.cc`, `prob_geom.cc`), perturbations, and output (`printing.cc`, `matlab_io.cc`)
- `src/noise/` — frequency-domain neutron noise solvers with complex arithmetic
- `src/rom/` — reduced-order models using POD/LUPOD snapshot methods

### Input Parameter Files (`.prm`)

All configuration is read via deal.II's `ParameterHandler` from `.prm` files. Key parameters:

| Parameter | Options |
|---|---|
| `Transport_Appr` | `Diffusion`, `SPN`, `SDPN`, `Full_SPN` |
| `Geometry_Type` | `Rectangular`, `Hexagonal`, `Unstructured`, `Composed` |
| `Solver_Type` | `bifpam`, `power_it`, `slepc_2g`, `slepc_7g`, `newton`, `gd`, `ks` |
| `Matrix_Free_Type` | `full_allocated`, `non_diagonal`, `full_matrixfree` |
| `Transient` | `true`/`false` — enables time-dependent mode |
| `Noise_Calculation` | `true`/`false` — enables frequency-domain noise |
| `ROM_Transient` / `ROM_Static` | `true`/`false` — enables ROM |

The `.prm` files in `examples/` subdirectories are the best reference for valid parameter combinations.

### Cross-Section Files

Materials are read from separate cross-section files (`.xs2g`, `.xsec`, `.xml`), referenced by `XSECS_Filename` in the `.prm`. The `Materials` class (`include/io/materials.h`) handles all XS I/O.

### Output

Results are written to text files and optionally `.vtk` files (for ParaView visualization). For noise calculations, `.mat` (MATLAB) output is supported via `RESULTS_Filename`.
