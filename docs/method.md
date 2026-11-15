# Discrete stationary Stokes model

The new solver is a **fresh implementation** grounded in the cavity boundary condition from `legacy/coursework/functioncall.m` and the vorticity–velocity–pressure setting of `legacy/coursework/project.m`. It intentionally does *not* assert equivalence to the original mixed finite-element code.

## Governing equations

On the unit square \(\Omega=(0,1)^2\), solve

\[
-\nu\Delta\mathbf{u}+\nabla p=\mathbf{f},\qquad
\nabla\cdot\mathbf{u}=0,
\]

with prescribed velocity on the boundary. \(\nu>0\) is nondimensional viscosity. There are **no nonlinear advection terms**. Therefore these are steady **Stokes equations**, not a finite-Reynolds-number Navier–Stokes simulation.

## MAC staggering

For \(n\) cells per side and \(h=1/n\):

- \(u_{i,j}\) at \((ih,(j-\tfrac12)h)\) for \(1\le i\le n-1\), \(1\le j\le n\): \(n(n-1)\) unknowns.
- \(v_{i,j}\) at \(((i-\tfrac12)h,jh)\) for \(1\le i\le n\), \(1\le j\le n-1\): \(n(n-1)\) unknowns.
- \(p_{i,j}\) at \(((i-\tfrac12)h,(j-\tfrac12)h)\) for \(1\le i,j\le n\): \(n^2\) unknowns, one pinned.

The momentum equations use a five-point discrete Laplacian; pressure gradients are differences of adjacent cell pressures divided by \(h\). Continuity is the signed sum of the four face-normal velocities divided by \(h\), with known boundary fluxes on the right-hand side.

Normal boundary velocities lie on actual MAC faces and are incorporated directly. Tangential velocities lie half a grid cell from walls and use mirror ghost values \(u_\mathrm{ghost}=2u_\mathrm{wall}-u_\mathrm{first}\). The same method applies to \(v\) at side walls.

## Block system and gauge

Let \(A\) be the sparse viscous operator and \(D\) the divergence of the **unknown** face velocities. With \(G=-D^T\) and prescribed boundary contribution \(b\), the equations are

\[
\begin{bmatrix}A&-D^T\\-D&0\end{bmatrix}
\begin{bmatrix}\mathbf{q}\\\mathbf{p}\end{bmatrix}
=
\begin{bmatrix}\mathbf{f}_\mathrm{adjusted}\\-\mathbf{b}\end{bmatrix}.
\]

Pressure is unique only up to a constant. We set the first cell pressure to zero, delete its column, and drop the corresponding redundant continuity row. We require the boundary's net normal flux to be zero. The solver reports the **full** continuity residual including that removed row.

The reduced system has \(2n(n-1)+n^2-1\) unknowns and is solved using MATLAB's sparse backslash. The implementation is meant for moderate grids, not optimized large-scale CFD.

## Verification cases

1. **Zero:** zero forcing and homogeneous velocity boundary values give zero velocity and pressure.
2. **Linear shear:** \(u=y, v=0,p=0\) with \(u=y\) on every boundary. This has balanced inflow/outflow at the side walls, and tests nonzero normal boundary conditions; the MAC discretization should reproduce it up to roundoff.
3. **Manufactured exact solution:**

   \[
   \psi(x,y)=\sin^2(\pi x)\sin^2(\pi y),\quad
   u=\partial_y\psi=\pi\sin^2(\pi x)\sin(2\pi y),
   \]

   \[
   v=-\partial_x\psi=-\pi\sin(2\pi x)\sin^2(\pi y),\quad
   p=\cos(\pi x)\cos(\pi y).
   \]

   The force is computed analytically as \(-\nu\Delta(u,v)+\nabla p\) in `stokes.make_case`. The pressure is compared after shifting to the same gauge.
4. **Lid-driven cavity:** \(u=1\) along the open top edge, and zero normal/tangential velocity at the other walls. The idealized corner discontinuity can generate pressure singular behavior. Only algebraic consistency, discrete incompressibility, zero penetration and nontrivial interior flow are asserted—not benchmark agreement against a finite-Reynolds-number problem.

## Important interpretation

- The velocity vectors displayed by `stokes.plot_solution` are **cell-center averages** of the staggered unknowns, used only for plotting.
- A tiny discrete divergence is an **algebraic conservation check**, not proof of physical accuracy. The manufactured-solution refinement test supplies an independent accuracy test.
- The code does not implement nonlinear convection, turbulence, Reynolds-number sweeps, or transient time integration.
