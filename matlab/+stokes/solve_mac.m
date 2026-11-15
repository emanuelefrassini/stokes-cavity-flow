function sol = solve_mac(n, problem)
%SOLVE_MAC Solve the stationary incompressible Stokes system (sparse direct).
%   sol = stokes.solve_mac(n, stokes.make_case('cavity',1));
sys = stokes.assemble_mac(n,problem);
z = sys.K\sys.rhs;
vel = z(1:sys.Nvel);
p = [0; z(sys.Nvel+1:end)];  % gauge: pressure at cell (1,1) is zero
h = sys.h;

u = reshape(vel(1:sys.Nu),n-1,n);
v = reshape(vel(sys.Nu+1:end),n,n-1);
p_center = reshape(p,n,n);

% Full boundary-normal face velocities (u: (n+1)x n, v: n x(n+1)).
ufaces = zeros(n+1,n);
vfaces = zeros(n,n+1);
ufaces(2:n,:) = u;
vfaces(:,2:n) = v;
for j = 1:n
    y = (j-0.5)*h;
    ufaces(1,j) = problem.u_boundary(0,y);
    ufaces(n+1,j) = problem.u_boundary(1,y);
end
for i = 1:n
    x = (i-0.5)*h;
    vfaces(i,1) = problem.v_boundary(x,0);
    vfaces(i,n+1) = problem.v_boundary(x,1);
end

% The following averages are for plots only; the actual unknowns are staggered.
uc = 0.5*(ufaces(1:n,:) + ufaces(2:n+1,:));
vc = 0.5*(vfaces(:,1:n) + vfaces(:,2:n+1));
[Xc,Yc] = ndgrid(((1:n)-0.5)*h,((1:n)-0.5)*h);
div = reshape(sys.D*vel - sys.b_div,n,n);
mom = sys.A*vel + sys.G*p - sys.f;
scale = max([1, norm(sys.f,inf), norm(sys.b_div,inf)]);

sol.n = n;
sol.h = h;
sol.nu = sys.nu;
sol.problem = problem.name;
sol.u = u;
sol.v = v;
sol.pressure = p_center;
sol.u_faces = ufaces;
sol.v_faces = vfaces;
sol.u_center = uc;
sol.v_center = vc;
sol.x_center = Xc;
sol.y_center = Yc;
sol.divergence = div;
sol.max_divergence = max(abs(div(:)));
sol.relative_residual = max([norm(mom,inf),norm(div(:),inf)])/scale;
sol.gauge_pressure = p(1);
sol.system = sys;
end
