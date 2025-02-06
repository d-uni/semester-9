% Parameters
h = 0.05; % Step size
T = 6;    % End time
N = T / h; % Number of steps
t = 0:h:T; % Time vector
y_exact_solution = @(t) exp(-5 * t); % Exact solution
f = @(t, y) -5 * y; % Define the function f(t, y)

% Initialize arrays for solutions
y_exact = y_exact_solution(t); % Exact solution at all points
y_simpson_exact = zeros(1, N + 1); % Solution using Simpson with exact y1
y_simpson_rk4 = zeros(1, N + 1);   % Solution using Simpson with RK4 y1
y_simpson_euler = zeros(1, N + 1); % Solution using Simpson with Euler y1

% Initial condition
y_simpson_exact(1) = 1;
y_simpson_rk4(1) = 1;
y_simpson_euler(1) = 1;

% Compute y1 using three methods
% (a) Exact solution
y1_exact = y_exact_solution(h);

% (b) 4th order Runge-Kutta method
k1 = f(t(1), y_simpson_rk4(1));
k2 = f(t(1) + h / 2, y_simpson_rk4(1) + h * k1 / 2);
k3 = f(t(1) + h / 2, y_simpson_rk4(1) + h * k2 / 2);
k4 = f(t(1) + h, y_simpson_rk4(1) + h * k3);
y1_rk4 = y_simpson_rk4(1) + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

% (c) Forward Euler method
y1_euler = y_simpson_euler(1) + h * f(t(1), y_simpson_euler(1));

% Assign y1 to solutions
y_simpson_exact(2) = y1_exact;
y_simpson_rk4(2) = y1_rk4;
y_simpson_euler(2) = y1_euler;

% Two-step Simpson's method
for n = 1:N-1
    %disp(y_simpson_exact(n))
    %disp(y_simpson_euler(n+1))

    %disp(y_simpson_euler(n))
    %disp(y_simpson_euler(n+1))
    %disp('-------------------')
    
    % (a) Using exact y1
    func_exact = @(y_next) y_simpson_exact(n) + (h / 3) * ...
        (f(t(n), y_simpson_exact(n)) + 4 * f(t(n+1), y_simpson_exact(n+1)) + f(t(n+2), y_next)) - y_next;
    % Solve for y_{n+2} using fsolve
    y_guess = y_simpson_exact(n+1); % Initial guess
    y_simpson_exact(n+2) = fsolve(func_exact, y_guess, optimset('Display', 'off'));
    
    % (b) Using RK4 y1
    func_rk4 = @(y_next) y_simpson_rk4(n) + (h / 3) * ...
        (f(t(n), y_simpson_rk4(n)) + 4 * f(t(n+1), y_simpson_rk4(n+1)) + f(t(n+2), y_next)) - y_next;
    y_guess = y_simpson_rk4(n+1); % Initial guess
    y_simpson_rk4(n+2) = fsolve(func_rk4, y_guess, optimset('Display', 'off'));
    
    % (c) Using Euler y1
    func_euler = @(y_next) y_simpson_euler(n) + (h / 3) * ...
        (f(t(n), y_simpson_euler(n)) + 4 * f(t(n+1), y_simpson_euler(n+1)) + f(t(n+2), y_next)) - y_next;
    y_guess = y_simpson_euler(n+1); % Initial guess
    y_simpson_euler(n+2) = fsolve(func_euler, y_guess, optimset('Display', 'off'));
end

% Compute errors
error_exact = abs(y_simpson_exact - y_exact);
error_rk4 = abs(y_simpson_rk4 - y_exact);
error_euler = abs(y_simpson_euler - y_exact);

% Plot the errors
figure;
semilogy(t, error_exact, 'b-', 'LineWidth', 1.5); hold on;
semilogy(t, error_rk4, 'r--', 'LineWidth', 1.5);
semilogy(t, error_euler, 'g-.', 'LineWidth', 1.5);

% Formatting the plot
xlabel('Time (t)');
ylabel('Error semilogy(|y_{numerical} - y_{exact}|)');
legend('Simpson (Exact y_1)', 'Simpson (RK4 y_1)', 'Simpson (Euler y_1)');
title('Error Comparison for Two-Step Simpsons Method');
grid on;

% Display final errors
fprintf('Final errors at T = %.2f:\n', T);
fprintf('Using exact y1: %.8f\n', error_exact(end));
fprintf('Using RK4 y1:   %.8f\n', error_rk4(end));
fprintf('Using Euler y1: %.8f\n', error_euler(end));
