%REFINEMENT_STUDY Verified manufactured Stokes solution on nested grids.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab'));
ns = [4,8,16,24,32];
problem = stokes.make_case('manufactured',1);
E = zeros(numel(ns),5);
for k = 1:numel(ns)
    sol = stokes.solve_mac(ns(k),problem);
    m = stokes.error_metrics(sol,problem);
    E(k,:) = [ns(k),m.velocity_max,m.pressure_max,m.divergence_max,m.linear_residual];
    fprintf('n=%2d velocity=%.4e pressure=%.4e div=%.3e residual=%.3e\n', E(k,:));
    if k > 1
        rv = log(E(k-1,2)/E(k,2))/log(ns(k)/ns(k-1));
        rp = log(E(k-1,3)/E(k,3))/log(ns(k)/ns(k-1));
        fprintf('   observed orders: velocity=%.3f, pressure=%.3f\n',rv,rp);
    end
end
outdir = fullfile(root,'results');
if ~exist(outdir,'dir')
    mkdir(outdir);
end
f = fopen(fullfile(outdir,'refinement.csv'),'w');
assert(f >= 0,'Cannot create results/refinement.csv');
fprintf(f,'n,velocity_max,pressure_max,divergence_max,relative_residual\n');
for k = 1:size(E,1)
    fprintf(f,'%d,%.16g,%.16g,%.16g,%.16g\n',E(k,:));
end
fclose(f);
