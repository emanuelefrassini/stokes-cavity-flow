`refinement_study.m` runs the exact smooth Stokes problem at five grid sizes and saves errors to `results/refinement.csv`.

For MATLAB, from the project root:

```matlab
run('experiments/refinement_study.m')
```

The script prints empirical convergence rates for velocity and pressure between successive grids.
