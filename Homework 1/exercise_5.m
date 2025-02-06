% Parameters for the Lotka-Volterra equations
alpha = 0.2;  % Growth rate of prey
beta = 0.01;  % Predation rate
gamma = 0.004; % Predator reproduction rate
delta = 0.07; % Predator mortality rate

% Initial conditions
x0 = 19;  % Initial prey population
y0 = 22;  % Initial predator population
t0 = 0;   % Start time
T = 300;  % End time
h = 1e-3; % Step size

% Number of steps
N = round((T - t0) / h);

% Initialize arrays for storing results
t = linspace(t0, T, N+1); % Time array
x = zeros(1, N+1);        % Prey population array
y = zeros(1, N+1);        % Predator population array

% Initial values
x(1) = x0;
y(1) = y0;

% Define the Lotka-Volterra system as functions
dx = @(x, y) x * (alpha - beta * y); % dx/dt
dy = @(x, y) y * (gamma * x - delta); % dy/dt

% Implement 4th-order Runge-Kutta method
for n = 1:N
    % Current values
    xn = x(n);
    yn = y(n);
    
    % RK4 coefficients for x
    k1x = h * dx(xn, yn);
    k2x = h * dx(xn + 0.5 * k1x, yn + 0.5 * h * dy(xn, yn));
    k3x = h * dx(xn + 0.5 * k2x, yn + 0.5 * h * dy(xn, yn));
    k4x = h * dx(xn + k3x, yn + h * dy(xn, yn));
    
    % RK4 coefficients for y
    k1y = h * dy(xn, yn);
    k2y = h * dy(xn + 0.5 * h * dx(xn, yn), yn + 0.5 * k1y);
    k3y = h * dy(xn + 0.5 * h * dx(xn, yn), yn + 0.5 * k2y);
    k4y = h * dy(xn + h * dx(xn, yn), yn + k3y);
    
    % Update populations
    x(n+1) = xn + (k1x + 2*k2x + 2*k3x + k4x) / 6;
    y(n+1) = yn + (k1y + 2*k2y + 2*k3y + k4y) / 6;
end

% Display final results
disp('Final populations:');
disp(['Prey: ', num2str(x(end))]);
disp(['Predators: ', num2str(y(end))]);

% Plot the results
figure;

% Plot prey population
plot(t, x, 'b', 'LineWidth', 1.5);
hold on;

% Plot predator population
plot(t, y, 'r', 'LineWidth', 1.5);

% Add labels, title, and legend
xlabel('Time');
ylabel('Population');
title('Lotka-Volterra Predator-Prey Model');
legend('Prey (x)', 'Predator (y)', 'Location', 'Best');
grid on;

% Show the plot
hold off;
