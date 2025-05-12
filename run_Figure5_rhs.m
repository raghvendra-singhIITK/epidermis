y0=[0.5,0.5,0.5, 0];
tspan=[0 6000];
[t, y]=ode45('Figure5_rhs', tspan, y0);

x1=y(:,2);
plot(t, x1, 'LineWidth', 2, 'color', 'k');

