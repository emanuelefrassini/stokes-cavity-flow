function problem = make_case(name, nu)
%MAKE_CASE Construct boundary data, forcing and optional exact solution.
%   problem = stokes.make_case('cavity', nu) produces a unit lid-speed
%   stationary Stokes problem on [0,1]^2. The lid is discontinuous at the
%   corners, as in the standard cavity idealization.
%
%   problem = stokes.make_case('manufactured', nu) gives a smooth,
%   divergence-free, no-slip exact solution with a known pressure.
%
%   No toolboxes are required. All coordinates are nondimensional.
if nargin < 2
    nu = 1.0;
end
if ~isscalar(nu) || ~isfinite(nu) || nu <= 0
    error('stokes:InvalidViscosity', 'nu must be a positive finite scalar.');
end
problem.nu = nu;
problem.name = lower(char(name));
problem.pressure_exact = [];

switch problem.name
    case {'cavity', 'lid', 'lid-driven-cavity'}
        problem.name = 'cavity';
        % The ghost-cell treatment samples the lid only at interior faces.
        problem.u_boundary = @(x,y) double(abs(y-1) < 1e-14) + 0*x;
        problem.v_boundary = @(x,y) 0*x + 0*y;
        problem.force_u = @(x,y) 0*x + 0*y;
        problem.force_v = @(x,y) 0*x + 0*y;
        problem.u_exact = [];
        problem.v_exact = [];
    case 'manufactured'
        % Streamfunction psi = sin(pi*x)^2 * sin(pi*y)^2.
        % u = d(psi)/dy, v = -d(psi)/dx, hence div(u) = 0.
        pi0 = pi;
        uf = @(x,y) pi0.*sin(pi0.*x).^2.*sin(2*pi0.*y);
        vf = @(x,y) -pi0.*sin(2*pi0.*x).*sin(pi0.*y).^2;
        pf = @(x,y) cos(pi0.*x).*cos(pi0.*y);
        problem.u_boundary = uf;
        problem.v_boundary = vf;
        problem.u_exact = uf;
        problem.v_exact = vf;
        problem.pressure_exact = pf;
        problem.force_u = @(x,y) -nu*2*pi0^3.*(2*cos(2*pi0.*x)-1).*sin(2*pi0.*y) ...
                                   -pi0.*sin(pi0.*x).*cos(pi0.*y);
        problem.force_v = @(x,y) -nu*2*pi0^3.*sin(2*pi0.*x).*(1-2*cos(2*pi0.*y)) ...
                                   -pi0.*cos(pi0.*x).*sin(pi0.*y);
    case {'linear-shear','shear'}
        % Exact u = y, v = 0, p = 0 (with balanced side throughflow).
        problem.name = 'linear-shear';
        problem.u_boundary = @(x,y) y + 0*x;
        problem.v_boundary = @(x,y) 0*x + 0*y;
        problem.force_u = @(x,y) 0*x + 0*y;
        problem.force_v = @(x,y) 0*x + 0*y;
        problem.u_exact = @(x,y) y + 0*x;
        problem.v_exact = @(x,y) 0*x + 0*y;
        problem.pressure_exact = @(x,y) 0*x + 0*y;
    case 'zero'
        problem.u_boundary = @(x,y) 0*x + 0*y;
        problem.v_boundary = @(x,y) 0*x + 0*y;
        problem.force_u = @(x,y) 0*x + 0*y;
        problem.force_v = @(x,y) 0*x + 0*y;
        problem.u_exact = @(x,y) 0*x + 0*y;
        problem.v_exact = @(x,y) 0*x + 0*y;
        problem.pressure_exact = @(x,y) 0*x + 0*y;
    otherwise
        error('stokes:UnknownCase', 'Unknown case: %s.', name);
end
end
