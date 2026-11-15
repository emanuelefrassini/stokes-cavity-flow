function fig = plot_solution(sol)
%PLOT_SOLUTION Display cell-centered interpolated velocity and pressure.
%   Cell-centered velocities are averages of staggered face values.
X = sol.x_center';
Y = sol.y_center';
U = sol.u_center';
V = sol.v_center';
P = sol.pressure';
S = hypot(U,V);
fig = figure('Color','w','Position',[100,100,1050,760]);
subplot(2,2,1)
contourf(X,Y,S,18,'LineStyle','none'); colorbar;
hold on;
stride = max(1,floor(sol.n/16));
quiver(X(1:stride:end,1:stride:end),Y(1:stride:end,1:stride:end), ...
       U(1:stride:end,1:stride:end),V(1:stride:end,1:stride:end),'k');
hold off; axis equal tight; xlabel('x'); ylabel('y'); title('Velocity magnitude and vectors');
subplot(2,2,2)
contourf(X,Y,U,18,'LineStyle','none'); colorbar;
axis equal tight; xlabel('x'); ylabel('y'); title('Horizontal velocity u');
subplot(2,2,3)
contourf(X,Y,P,18,'LineStyle','none'); colorbar;
axis equal tight; xlabel('x'); ylabel('y'); title('Pressure (p_{1,1} = 0)');
subplot(2,2,4)
contourf(X,Y,sol.divergence',18,'LineStyle','none'); colorbar;
axis equal tight; xlabel('x'); ylabel('y'); title('Cell-wise continuity residual');
end
