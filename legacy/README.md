# Original MATLAB coursework (historical, unverified)

These scripts are preserved **unaltered** from the uploaded CFD coursework, for provenance and comparison; they are **not** required to run the new solver.

- `project.m` / `functioncall.m`: original vorticity–velocity–pressure mixed-system experiment, including a lid-driven boundary condition.
- `new_project_test.m`: an alternative draft; its first declaration names `project` although the filename is different, and it overrides its grid argument internally. It should not be on the active MATLAB path.
- `simple.m`, `PROJECT2.m`, `project2try.m`, `project2_1.m`: related pressure–velocity and block-system experiments. They are not equivalent to the new MAC discretization and are not validated here.
- `gauss_seidel.m`: historical iterative experiment, not called by the verified solver.

## Issues noted in the historic `project.m`

1. The velocity and pressure arrays are resized for visualization using sums of adjacent values rather than a documented staggered interpolation; these should not be interpreted as solved boundary values.
2. The expression `phi = v1try.*y-v2try.*x` is **not generally a streamfunction**; a streamfunction must satisfy `u = dpsi/dy` and `v = -dpsi/dx`.
3. The original matrix construction uses dense `diag`, `kron`, and `zeros` blocks and deletes columns/rows repeatedly. It becomes expensive on refined meshes.
4. The original assembled system has not been independently checked for consistency, pressure uniqueness, or grid convergence. The current project **does not claim** that those scripts produced a validated Stokes solution.

To avoid function-name collisions, do **not** add this directory to the MATLAB path or call `addpath(genpath(pwd))`. The project runners add `matlab/` explicitly.

Only user-supplied coursework files were copied. Third-party source files from the larger CFD archive were not included.
