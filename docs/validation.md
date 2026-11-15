# Numerical validation and traceability

## GNU Octave execution (user-reported, 15 November 2026)

`run_tests` passed under GNU Octave on WSL. The successful cases in the suite cover invalid-parameter handling, sparse matrix assembly, pressure gauge and adjoint structure, zero flow, exact linear shear, lid-driven cavity continuity, and convergence to a smooth manufactured solution. These are **numerical regression tests**, not a full CFD benchmark against physical measurements.

`run_demo` reported for the steady Stokes cavity with `n=24`, `nu=1`:

```text
Steady Stokes cavity: n=24, nu=1
Unknowns: 1679 (velocity 1104 + pressure 575 minus gauge)
Maximum discrete divergence: 4.9735e-15
Relative saddle-system residual: 3.9475e-16
```

The following are exact transcriptions of the user-reported grid refinement summary (rounded to the displayed digits):

| n | Max velocity error | Max pressure error | Max discrete divergence | Relative system residual | Observed velocity order | Observed pressure order |
|---:|---:|---:|---:|---:|---:|---:|
| 4  | 5.1915e-01 | 4.4679e-02 | 1.776e-15 | 2.114e-16 | — | — |
| 8  | 1.5392e-01 | 1.2418e-02 | 3.553e-15 | 6.515e-16 | 1.754 | 1.847 |
| 16 | 3.9904e-02 | 3.1855e-03 | 1.188e-14 | 2.921e-15 | 1.948 | 1.963 |
| 24 | 1.7851e-02 | 1.4225e-03 | 2.842e-14 | 8.368e-15 | 1.984 | 1.988 |
| 32 | 1.0064e-02 | 8.0148e-04 | 8.711e-14 | 1.046e-14 | 1.992 | 1.994 |

Errors and estimated orders correspond to the manufactured solution described in [method.md](method.md). The result is consistent with approximately second-order spatial convergence. `experiments/refinement_study.m` generates a new, higher-precision `results/refinement.csv` on the machine where it is run; the CSV in `docs/` stores only the rounded console observations shown above.

The `run_demo` command also emitted libEGL/Mesa messages while attempting to initialize accelerated graphics under WSL. The **numerical solve completed**. The existence and quality of the saved PNG should be checked separately, for example with `ls -lh figures/cavity_n24.png` and by opening the image. Graphics warnings should not be interpreted as a failed Stokes solve.

GNU Octave version was not supplied with the log. MATLAB itself has **not** been executed in the current review.

## Independent mathematical cross-check

Before the MATLAB-compatible implementation was prepared, an independently written Python/SciPy sparse discretization of the same MAC equations was used to check the intended matrix assembly and grid convergence. Those results were a design cross-check, not evidence that MATLAB files had executed. The subsequent GNU Octave run now supplies implementation-level numerical verification.

## CI validation

`.github/workflows/matlab-tests.yml` installs **GNU Octave** on an Ubuntu runner and executes `run_tests`. No MATLAB licence is required for CI. The CI run still must complete successfully on GitHub after publication; local tests do not imply CI is green.

## Scientific limitations

The Stokes cavity calculation has no nonlinear convective term; it is not a Navier–Stokes benchmark for finite Reynolds number. Near-zero algebraic residuals demonstrate consistency with the **discrete** model, not experimental accuracy. The original mixed-method coursework scripts in `legacy/coursework/` remain historical and unverified.
