%a) --------------------------------------------
% Define the stability function R(z)
R = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24; % Stability function R(z)

% Define the function to solve |R(z)| = 1
stability_condition = @(z) abs(R(z)) - 1; % |R(z)| = 1 condition

% Use fminsearch to find the root where |R(z)| = 1
z_initial_guess = -2; % Initial guess for the search
z_solution = fminsearch(@(z) abs(stability_condition(z)), z_initial_guess);

% Display the solution
disp(['The boundary of the stability region is approximately z = ', num2str(z_solution)]);


% Define problem parameters
nx = 60;
G = numgrid('S', nx);
A = delsq(G) * (nx - 1)^2; % Generate A matrix
y0 = ones(size(A, 1), 1);  % Initial condition
T = 0.1;                  % Final time
h_values = [1e-2, 1e-3, 1e-4]; % Step sizes

%Compute eigenvalue with largest modulus
lambda = -eigs(A, 1, 'lm');

% Compute interval of absolute stability
% For the 4th-order Runge-Kutta method, the stability region is approximately:
% |lambda * h| <= -2.7853 (RK4 stability limit)
h_max = z_solution / abs(lambda);

disp(['Maximum step size for absolute stability: h <= ', num2str(h_max)]);

% Solve the ODE numerically using ode45
tic; % Start timing
[t, y_approx] = ode45(@(t, y) -A * y, [0 T], y0);
cpu_time = toc; % End timing

% Compute exact solution
y_exact = expm(-T * A) * y0;

% Step 5: Compute the infinity norm of the error
error_norm = norm(y_exact - y_approx(end, :)', inf);

% Display results
disp(['Infinity norm of the error: ', num2str(error_norm)]);
disp(['CPU time elapsed: ', num2str(cpu_time), ' seconds']);



%b)-------------------------------------------

% Initialize results table
results = [];

% Crank-Nicolson Method
for h = h_values
    tic; % Start timing
    N = round(T / h); % Number of time steps
    y_CN = y0;
    I = speye(size(A));
    for k = 1:N
        % Solve (I + h/2 * A) y_(n+1) = (I - h/2 * A) y_n
        b = (I - h/2 * A) * y_CN;
        y_CN = cgs((I + h/2 * A), b, h^3); % Solve system of linear equations with \
    end
    cpu_time = toc; % End timing
    error_norm = norm(y_exact - y_CN, inf); % Infinity norm of the error
    results = [results; {'Crank-Nicolson', h, N, error_norm, cpu_time}];
end

% BDF3 Method
for h = h_values
    tic; % Start timing
    N = round(T / h); % Number of time steps
    I = speye(size(A)); % Identity matrix of the same size as A
    % Start with first two steps using Crank Nicolson method
    y_old = y0; % Initial y(0)
    y_BDF3 = (I + h/2 * A) \ ((I - h/2 * A) * y0); % Second step with CN
    y_old2 = y_BDF3;
    y_BDF3 = (I + h/2 * A) \ ((I - h/2 * A) * y_old2); % Solve for second step
    % Main BDF3 loop
    for k = 3:N
        % Set up the right-hand side of the BDF3 equation
        b = (3 * y_BDF3 - 3/2 * y_old2 + 1/3 * y_old);
        % Update old values
        y_old = y_old2;
        y_old2 = y_BDF3;
        
        % Solve the BDF3 equation for y_(n+1)
        y_BDF3 = cgs((11/6 * I + h * A), b, h^3);
    end

    % End timing
    cpu_time = toc;

    % Calculate the error (difference from exact solution at final time)
    error_norm = norm(y_exact - y_BDF3, inf); % Infinity norm of the error

    % Store results
    results = [results; {'BDF3', h, N, error_norm, cpu_time}];
end

% Display results table
disp('Results Table:');
disp(table(results(:, 1), cell2mat(results(:, 2)), cell2mat(results(:, 3)), ...
           cell2mat(results(:, 4)), cell2mat(results(:, 5)), ...
           'VariableNames', {'Method', 'StepSize', 'Steps', 'Error', 'CPU_Time'}));