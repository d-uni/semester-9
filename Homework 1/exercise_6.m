clc;
clear;

% Parameters
y0 = 1;            % Initial condition
T = 100;           % Final time
h_init = 1e-3;     % Initial timestep for adaptive method
h_const = 1e-3;    % Constant timestep
tol = 1e-8;        % Tolerance for adaptive method
f = @(y) -10*(y^2);

% Exact solution
exact_solution = @(t) 1 ./ (1 + 10 * t);

% Crank-Nicolson step function
crank_nicolson = @(y_prev, h) fzero(@(y_next) ...
    y_next - y_prev - (h/2) * (f(y_prev) + f(y_next)), y_prev);

% Constant Timestep Solution
t_const = 0:h_const:T;
y_const = zeros(size(t_const));
y_const(1) = y0;

for i = 2:length(t_const)
    y_const(i) = crank_nicolson(y_const(i-1), h_const);
end


% Compute error for constant timestep solution
exact_const = exact_solution(t_const);
error_const = abs(y_const - exact_const);

% Adaptive Timestep Solution
t_adaptive = 0;
y_adaptive = y0;
h_adaptive = h_init;
h_values = h_adaptive;  % Store step sizes
errors_adaptive = [];   % Store errors for adaptive method

while t_adaptive(end) + h_adaptive < T

    % Full-step solution
    y_full = crank_nicolson(y_adaptive(end), h_adaptive);

    % Half-step solution
    y_half_1 = crank_nicolson(y_adaptive(end), h_adaptive / 2);
    y_half_2 = crank_nicolson(y_half_1, h_adaptive / 2);

    % Error estimate
    error_estimate = (4 / 3) * abs(y_full - y_half_2);

    % Adjust timestep
    h_adaptive = h_adaptive * ((tol / error_estimate)^(1/3));

    t_adaptive(end + 1) = t_adaptive(end) + h_adaptive;
    y_adaptive(end + 1) = crank_nicolson(y_adaptive(end), h_adaptive);
    errors_adaptive(end + 1) = error_estimate;
    h_values(end + 1) = h_adaptive;    
        
    
end
h_adaptive = T - t_adaptive(end);
y_full = crank_nicolson(y_adaptive(end), h_adaptive);
t_adaptive(end + 1) = t_adaptive(end) + h_adaptive;
y_adaptive(end + 1) = y_full;
h_values(end + 1) = h_adaptive;

% Compute exact solution and error for adaptive method
exact_adaptive = exact_solution(t_adaptive);
error_adaptive = abs(y_adaptive - exact_adaptive);

% Plots
% 1. Error comparison plot
figure;
semilogy(t_const, error_const, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Constant Timestep');
hold on;
semilogy(t_adaptive, error_adaptive, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Adaptive Timestep');
xlabel('Time');
ylabel('Error');
title('Error with Reference to Exact Solution');
legend('Location', 'northeast');
grid on;

% 2. Timestep size plot
figure;
plot(t_adaptive(2:end), h_values(2:end), 'k-', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Timestep Size (h)');
title('Timestep Size vs. Time (Adaptive Method)');
grid on;

