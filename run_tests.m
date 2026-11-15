function run_tests()
%RUN_TESTS Toolbox-free numerical regression tests for MAC Stokes solver.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'matlab'));

% Grid/parameter validation.
failed = false;
try
    stokes.make_case('cavity',0);
catch ME
    failed = strcmp(ME.identifier,'stokes:InvalidViscosity');
end
assert(failed,'Expected invalid-viscosity error.');
failed = false;
try
    stokes.assemble_mac(2,stokes.make_case('zero'));
catch ME
    failed = strcmp(ME.identifier,'stokes:InvalidGrid');
end
assert(failed,'Expected invalid-grid error.');
fprintf('PASS: parameter validation\n');

% Sparse block dimensions, adjoints and positive velocity block.
sys = stokes.assemble_mac(8,stokes.make_case('cavity'));
assert(issparse(sys.K) && issparse(sys.A) && issparse(sys.D));
assert(isequal(size(sys.A),[112,112]));
assert(isequal(size(sys.D),[64,112]));
assert(isequal(size(sys.K),[175,175]));
assert(norm(sys.G+sys.D','fro') < 1e-13);
assert(norm(sys.K-sys.K','fro') < 1e-13);
assert(all(diag(chol(sys.A))>0));
% Unit lid speed is imposed through tangential ghost-cell contributions.
n = 8;
first_top_u = (n-1)*(n-1)+1;
last_top_u = n*(n-1);
assert(max(abs(sys.f(first_top_u:last_top_u)-2*n*n)) < 1e-12);
fprintf('PASS: dimensions, symmetry, adjoints, positive velocity block, lid forcing\n');

% Homogeneous problem must give exactly zero solution (up to roundoff).
z = stokes.solve_mac(8,stokes.make_case('zero'));
assert(max(abs(z.u(:))) < 1e-13);
assert(max(abs(z.v(:))) < 1e-13);
assert(max(abs(z.pressure(:))) < 1e-13);
assert(z.max_divergence < 1e-13);
fprintf('PASS: zero-forcing zero-flow reference\n');

% Linear-shear manufactured solution is reproduced to roundoff.
linear_case = stokes.make_case('linear-shear',1);
for n = [4,8,12]
    s = stokes.solve_mac(n,linear_case);
    err = stokes.error_metrics(s,linear_case);
    assert(err.velocity_max < 1e-10);
    assert(err.pressure_max < 1e-10);
    assert(err.divergence_max < 1e-10);
end
fprintf('PASS: exact linear shear and balanced boundary throughflow\n');

% Cavity: no penetration, no body force, nonzero interior motion.
for n = [8,16]
    sol = stokes.solve_mac(n,stokes.make_case('cavity'));
    assert(sol.max_divergence < 1e-9);
    assert(sol.relative_residual < 1e-9);
    assert(all(abs(sol.u_faces(1,:)) < 1e-12));
    assert(all(abs(sol.u_faces(end,:)) < 1e-12));
    assert(all(abs(sol.v_faces(:,1)) < 1e-12));
    assert(all(abs(sol.v_faces(:,end)) < 1e-12));
    assert(max(abs(sol.u(:))) > 0.05);
    assert(abs(sol.gauge_pressure) < 1e-12);
end
fprintf('PASS: cavity continuity, gauge, nontrivial flow, no penetration\n');

% Smooth manufactured solution: observed second-order velocity/pressure.
ns = [8,16,32];
Ev = zeros(size(ns));
Ep = zeros(size(ns));
problem = stokes.make_case('manufactured');
for k = 1:numel(ns)
    sol = stokes.solve_mac(ns(k),problem);
    metrics = stokes.error_metrics(sol,problem);
    Ev(k) = metrics.velocity_max;
    Ep(k) = metrics.pressure_max;
    assert(metrics.linear_residual < 1e-9);
    assert(metrics.divergence_max < 1e-9);
end
order_v = log(Ev(1:end-1)./Ev(2:end))/log(2);
order_p = log(Ep(1:end-1)./Ep(2:end))/log(2);
assert(all(order_v > 1.7), 'Velocity convergence below expected order.');
assert(all(order_p > 1.7), 'Pressure convergence below expected order.');
assert(Ev(end) < 0.011 && Ep(end) < 0.001);
fprintf('PASS: manufactured velocity and pressure (orders: %.3f, %.3f)\n', ...
    min(order_v),min(order_p));
fprintf('All MAC Stokes numerical tests PASSED.\n');
end
