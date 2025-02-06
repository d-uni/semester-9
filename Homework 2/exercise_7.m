clc;
clear;

A = load('mat13041.rig'); 
A = spconvert(A);

n = size(A, 1);
x_exact = 1 ./ sqrt((1:n)');

b = A * x_exact;
setup.type = 'crout';        
setup.droptol = 0.01;     
[L, U] = ilu(A, setup);

% Parameters for GMRES
restarts = [10, 20, 30, 50]; % Values of restart

% Initialize arrays to store results
num_iterations = zeros(length(restarts), 1);
relative_residuals = zeros(length(restarts), 1);
cpu_times = zeros(length(restarts), 1);

% Convergence profiles for all restarts
figure;
hold on;
colors = ['r', 'g', 'b', 'm']; % Colors for each restart value

tol = 1e-12;       % Convergence tolerance
maxit = 550;       % Maximum number of iterations
x0 = zeros(n, 1);  % Initial guess

for i = 1:length(restarts)
    restart = restarts(i);

    % Timer start
    tic;
    
    % Anonymous function to use preconditioner in GMRES
    [x_approx, flag, relres, iter, resvec] = gmres(A, b, restart, tol, maxit,L, U);
    % Timer end 
    elapsed_time = toc;
   
    % Compute total number of iterations
    total_iterations = (iter(1) - 1) * restart + iter(2);

    % Store results
    num_iterations(i) = total_iterations; % Total iterations for GMRES (sum over inner iterations)
    relative_residuals(i) = relres; % Relative residual at convergence
    cpu_times(i) = elapsed_time; % Time in seconds

    % Plot convergence profile
    plot(0:total_iterations, log(resvec), 'Color', colors(i), 'DisplayName', sprintf('Restart = %d', restart));
end

% Configure plot
xlabel('Iteration number');
ylabel('Relative residual');
title('Convergence profiles for GMRES with ILU(0.01) preconditioner');
legend('show');
grid on;
hold off;

% Display results in a table format
results_table = table(restarts', num_iterations, relative_residuals, cpu_times, ...
    'VariableNames', {'Restart', 'Total_Iterations', 'Relative_Residual', 'CPU_Time'});

% Print results to the command window
disp('Results for GMRES with ILU(0.01) preconditioner:');
disp(results_table);
