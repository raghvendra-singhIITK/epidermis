y0=[0.5,0.5,0.5, 9.613];
tspan=[0 8000];
[t, y]=ode45('Figure3_rhs', tspan, y0);
x1=y(:,1);

plot(t, x1, 'LineWidth', 2, 'color', 'k');