function sys = assemble_mac(n, problem)
%ASSEMBLE_MAC Assemble steady Stokes equations on a staggered MAC grid.
%   -nu*Delta(u) + grad(p) = f,  div(u) = 0 on (0,1)^2.
%   u is stored at vertical faces ((i*h),(j-1/2)*h),
%   v at horizontal faces (((i-1/2)*h),j*h), p at cell centers.
%   Dirichlet normal velocities are inserted into divergence RHS.
%   Tangential velocities use second-order ghost-cell elimination.
%
%   The pressure has a constant null mode. The solver removes the first
%   pressure degree of freedom and its redundant continuity equation.
if ~isscalar(n) || ~isfinite(n) || n ~= floor(n) || n < 3
    error('stokes:InvalidGrid', 'n must be an integer >= 3.');
end
nu = problem.nu;
h = 1/n;
c = nu/h^2;
Nu = (n-1)*n;
Nvel = 2*Nu;
Np = n*n;

% Each velocity row has at most five matrix nonzeros.
ia = zeros(5*Nvel,1);
ja = zeros(5*Nvel,1);
va = zeros(5*Nvel,1);
na = 0;
f = zeros(Nvel,1);

% Horizontal velocity at vertical grid faces.
for j = 1:n
    y = (j-0.5)*h;
    for i = 1:n-1
        x = i*h;
        k = u_index(i,j,n);
        f(k) = problem.force_u(x,y);
        diagval = 4*c;
        if i > 1
            putA(k,u_index(i-1,j,n),-c);
        else
            f(k) = f(k) + c*problem.u_boundary(0,y);
        end
        if i < n-1
            putA(k,u_index(i+1,j,n),-c);
        else
            f(k) = f(k) + c*problem.u_boundary(1,y);
        end
        if j > 1
            putA(k,u_index(i,j-1,n),-c);
        else
            diagval = diagval + c;
            f(k) = f(k) + 2*c*problem.u_boundary(x,0);
        end
        if j < n
            putA(k,u_index(i,j+1,n),-c);
        else
            diagval = diagval + c;
            f(k) = f(k) + 2*c*problem.u_boundary(x,1);
        end
        putA(k,k,diagval);
    end
end

% Vertical velocity at horizontal grid faces.
for j = 1:n-1
    y = j*h;
    for i = 1:n
        x = (i-0.5)*h;
        k = v_index(i,j,n,Nu);
        f(k) = problem.force_v(x,y);
        diagval = 4*c;
        if i > 1
            putA(k,v_index(i-1,j,n,Nu),-c);
        else
            diagval = diagval + c;
            f(k) = f(k) + 2*c*problem.v_boundary(0,y);
        end
        if i < n
            putA(k,v_index(i+1,j,n,Nu),-c);
        else
            diagval = diagval + c;
            f(k) = f(k) + 2*c*problem.v_boundary(1,y);
        end
        if j > 1
            putA(k,v_index(i,j-1,n,Nu),-c);
        else
            f(k) = f(k) + c*problem.v_boundary(x,0);
        end
        if j < n-1
            putA(k,v_index(i,j+1,n,Nu),-c);
        else
            f(k) = f(k) + c*problem.v_boundary(x,1);
        end
        putA(k,k,diagval);
    end
end
A = sparse(ia(1:na),ja(1:na),va(1:na),Nvel,Nvel);

% D maps unknown face velocities to cell-centered divergence.
% b_div accounts for the prescribed normal velocities on the outer walls.
id = zeros(4*Np,1);
jd = zeros(4*Np,1);
vd = zeros(4*Np,1);
nd = 0;
b_div = zeros(Np,1);
for j = 1:n
    y = (j-0.5)*h;
    for i = 1:n
        x = (i-0.5)*h;
        k = p_index(i,j,n);
        if i > 1
            putD(k,u_index(i-1,j,n),-1/h);
        else
            b_div(k) = b_div(k) + problem.u_boundary(0,y)/h;
        end
        if i < n
            putD(k,u_index(i,j,n),1/h);
        else
            b_div(k) = b_div(k) - problem.u_boundary(1,y)/h;
        end
        if j > 1
            putD(k,v_index(i,j-1,n,Nu),-1/h);
        else
            b_div(k) = b_div(k) + problem.v_boundary(x,0)/h;
        end
        if j < n
            putD(k,v_index(i,j,n,Nu),1/h);
        else
            b_div(k) = b_div(k) - problem.v_boundary(x,1)/h;
        end
    end
end
D = sparse(id(1:nd),jd(1:nd),vd(1:nd),Np,Nvel);

if abs(sum(b_div)) > 1e-12*max(1,norm(b_div,1))
    error('stokes:IncompatibleBoundaryFlux', ...
        'Net prescribed boundary flux must be zero for incompressible flow.');
end
% Grad(p) = -D' p, giving symmetric saddle-point blocks [A -D'; -D 0].
G = -D';
K = [A, G(:,2:end); -D(2:end,:), sparse(Np-1,Np-1)];
rhs = [f; -b_div(2:end)];

sys.n = n;
sys.h = h;
sys.nu = nu;
sys.Nu = Nu;
sys.Nvel = Nvel;
sys.Np = Np;
sys.A = A;
sys.D = D;
sys.G = G;
sys.f = f;
sys.b_div = b_div;
sys.K = K;
sys.rhs = rhs;

    function putA(row,col,value)
        na = na+1;
        ia(na) = row;
        ja(na) = col;
        va(na) = value;
    end
    function putD(row,col,value)
        nd = nd+1;
        id(nd) = row;
        jd(nd) = col;
        vd(nd) = value;
    end
end

function k = u_index(i,j,n)
k = i + (j-1)*(n-1);
end
function k = v_index(i,j,n,Nu)
k = Nu + i + (j-1)*n;
end
function k = p_index(i,j,n)
k = i + (j-1)*n;
end
