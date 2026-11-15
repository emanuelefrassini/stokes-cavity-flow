%RUN_DEMO Solve and visualize lid-driven cavity Stokes flow.
% Run from the repository root in MATLAB: run_demo
addpath(fullfile(fileparts(mfilename('fullpath')),'matlab'));
n = 24;
nu = 1.0;
problem = stokes.make_case('cavity',nu);
sol = stokes.solve_mac(n,problem);
fprintf('Steady Stokes cavity: n=%d, nu=%g\n',n,nu);
fprintf('Unknowns: %d (velocity %d + pressure %d minus gauge)\n', ...
    size(sol.system.K,1),sol.system.Nvel,sol.system.Np-1);
fprintf('Maximum discrete divergence: %.4e\n',sol.max_divergence);
fprintf('Relative saddle-system residual: %.4e\n',sol.relative_residual);
fig = stokes.plot_solution(sol);
if ~exist('figures','dir')
    mkdir('figures');
end
saveas(fig,fullfile('figures','cavity_n24.png'));
