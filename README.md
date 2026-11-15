# Stokes Flow and Lid-Driven Cavity — MATLAB / GNU Octave

**A toolbox-free, sparse finite-difference solver for the two-dimensional steady incompressible Stokes equations, with numerical verification and a lid-driven cavity example.**

I revisited my Computational Fluid Dynamics coursework and developed this standalone MATLAB-compatible project to explore staggered-grid discretizations, incompressibility constraints, pressure gauge fixing, and numerical convergence. I have preserved selected original coursework scripts under [`legacy/coursework/`](legacy/README.md), separately from the new, tested MAC-grid implementation. **The original mixed vorticity–velocity–pressure coursework solver has not been verified or claimed equivalent to this implementation.**

The code also runs in **free GNU Octave**, so a MATLAB licence is not required.

## Mathematical problem

On the unit square, solve the **steady Stokes equations**

\[
-\nu\Delta\mathbf{u}+\nabla p=\mathbf{f},
\qquad \nabla\cdot\mathbf{u}=0.
\]

The lid-driven cavity uses a unit horizontal velocity on the top wall and stationary other walls. This is a **creeping-flow/Stokes** problem, **not** the nonlinear finite-Reynolds-number Navier–Stokes cavity benchmark.

## Numerical method

- Uniform staggered **MAC grid**: face-centred velocities, cell-centred pressures.
- Sparse five-point viscous operator and conservative discrete divergence.
- Dirichlet velocity conditions, with mirrored ghost values for tangential components.
- One pressure degree of freedom pinned to eliminate the pressure null space.
- Sparse linear solution using MATLAB/Octave backslash.
- Verification against exact zero flow, exact linear shear, and a smooth manufactured solution.
- Visualization of velocity magnitude, vectors, pressure, and the discrete continuity residual.

See the [mathematical formulation](docs/method.md) for grid indexing, block matrices, and details of the boundary treatment.

## Quick start — GNU Octave (no MATLAB licence required)

From the repository root, after installing GNU Octave:

```bash
octave --no-gui --quiet --eval "run_tests"
octave --no-gui --quiet --eval "run_demo"
octave --no-gui --quiet --eval "run('experiments/refinement_study.m')"
```

On Ubuntu / WSL, GNU Octave can be installed with `sudo apt update && sudo apt install octave`.

In MATLAB, open the repository as your current directory and run:

```matlab
run_tests
run_demo
run('experiments/refinement_study.m')
```

`run_demo` saves `figures/cavity_n24.png`; the refinement study writes `results/refinement.csv`. Run-generated plots and intermediate numerical files are gitignored until intentionally reviewed and added. Some WSL configurations emit Mesa/EGL graphics warnings; inspect whether the PNG was saved before diagnosing a plotting failure.

## Numerical validation

The new MATLAB-compatible solver was **executed successfully in GNU Octave on WSL** on 15 November 2026. The user-reported `run_tests` suite passed. For the 24-by-24 cavity case, `run_demo` reported:

| Diagnostic | Value |
|---|---:|
| Unknowns | 1,679 |
| Maximum discrete divergence | 4.9735 × 10⁻¹⁵ |
| Relative saddle-system residual | 3.9475 × 10⁻¹⁶ |

The smooth manufactured-solution refinement study produced:

| Cells per side | Max velocity error | Max pressure error | Velocity order | Pressure order |
|---:|---:|---:|---:|---:|
| 4 | 5.1915e-1 | 4.4679e-2 | — | — |
| 8 | 1.5392e-1 | 1.2418e-2 | 1.754 | 1.847 |
| 16 | 3.9904e-2 | 3.1855e-3 | 1.948 | 1.963 |
| 24 | 1.7851e-2 | 1.4225e-3 | 1.984 | 1.988 |
| 32 | 1.0064e-2 | 8.0148e-4 | **1.992** | **1.994** |

The observed convergence approaches second order for both velocity and pressure. A small algebraic residual verifies the **discrete** system and mass conservation; the separate exact-solution study supports solution accuracy. These results are *user-provided GNU Octave output*, not a MATLAB run or a finite-Reynolds-number CFD benchmark. See [validation and scope](docs/validation.md) and [machine-readable observations](docs/octave_refinement_2026-11-15.csv).

## Repository structure

```text
matlab/+stokes/       Sparse MAC assembly, solver, cases, metrics, plots
run_tests.m           Numerical regression tests (toolbox-free)
run_demo.m            Steady Stokes cavity demonstration
experiments/          Manufactured-solution grid refinement
docs/                 Mathematics, verification and recorded Octave outputs
legacy/coursework/    Preserved original unverified MATLAB drafts
.github/workflows/    GNU Octave CI running numerical tests
```

GitHub Actions runs `run_tests` **with GNU Octave**, without requiring MATLAB credentials or a licence. CI tests the numerical solver, not the graphics display. The historical coursework sources are not on the active MATLAB/Octave search path.

## Scope and limitations

- **Steady Stokes only:** no transient flow, convective Navier–Stokes terms, turbulence model, or Reynolds-number sweeps.
- The unit-speed lid has corner discontinuities. Reported cavity checks test algebraic consistency and continuity, not a comparison to experimental or classic finite-Reynolds-number cavity profiles.
- The unmodified historical coursework scripts are included for provenance and have **not** been validated against the new solver.
- The plotted cell-centre velocity is interpolated from the staggered face velocities.
- MATLAB execution and exact version equivalence with GNU Octave have not been checked independently; the GNU Octave CLI numerical tests did pass.

## Attribution

This project originates in my CFD coursework. I retained selected original MATLAB scripts as historical references and organized a new sparse MAC-grid implementation with independent verification cases. Third-party coursework utilities are not redistributed. If publishing the repository under an open-source licence, choose one that matches your intended reuse permissions.
