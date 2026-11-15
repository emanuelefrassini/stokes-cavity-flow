function metrics = error_metrics(sol, problem)
%ERROR_METRICS Staggered-grid max errors for the manufactured solution.
if isempty(problem.u_exact) || isempty(problem.v_exact) || ...
        isempty(problem.pressure_exact)
    error('stokes:NoExactSolution','This problem has no supplied exact solution.');
end
n = sol.n;
h = sol.h;
[Xu,Yu] = ndgrid((1:n-1)*h,((1:n)-0.5)*h);
[Xv,Yv] = ndgrid(((1:n)-0.5)*h,(1:n-1)*h);
ue = problem.u_exact(Xu,Yu);
ve = problem.v_exact(Xv,Yv);
pe = problem.pressure_exact(sol.x_center,sol.y_center);
% Compare pressures in the same gauge (p(1,1) = 0).
pe = pe-pe(1,1);
metrics.u_max = max(abs(sol.u(:)-ue(:)));
metrics.v_max = max(abs(sol.v(:)-ve(:)));
metrics.velocity_max = max(metrics.u_max,metrics.v_max);
metrics.pressure_max = max(abs(sol.pressure(:)-pe(:)));
metrics.divergence_max = sol.max_divergence;
metrics.linear_residual = sol.relative_residual;
end
