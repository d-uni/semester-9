% Parameters
T = 2; % Final time
k_values = 5:10; % Values of k for step size calculation
h_values = 2.^(-k_values); % Step sizes
exact_solution = @(t) 1 ./ (1 + 10 * t); % Exact solution function

% Preallocate for results
errors = zeros(size(h_values));

% Define the ODE function
f = @(t, y) -10 * y^2;

% 4-stage Runge-Kutta method
function y_next = runge_kutta_4_step(f, t, y, h)
    k1 = f(t, y);
    k2 = f(t + h/2, y + h/2 * k1);
    k3 = f(t + h/2, y + h/2 * k2);
    k4 = f(t + h, y + h * k3);
    y_next = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% Loop over different step sizes
for i = 1:length(h_values)
    h = h_values(i);
    N = T / h; % Number of steps
    t = 0:h:T; % Time vector
    y = 1; % Initial condition
    
    % Perform RK4 integration
    for n = 1:N
        y = runge_kutta_4_step(f, t(n), y, h);
    end
    
    % Compute error at T
    y_exact = exact_solution(T);
    errors(i) = abs(y - y_exact);
end

% Create log-log plot of errors
figure;
loglog(T./h_values, errors, 'o-', 'LineWidth', 1.5);
xlabel('Number of Steps');
ylabel('Error at T = 2');
title('Error vs Step Size (h) for RK4 Method');
grid on;

% Display results in table form
results_table = table(h_values', errors', 'VariableNames', {'Step_Size_h', 'Error_at_T'});
disp(results_table);
