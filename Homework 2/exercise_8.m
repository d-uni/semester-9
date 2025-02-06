clc;
clear;

% Load and convert the matrix from .mtx to MATLAB sparse format
A = load('ML_laplace.mtx'); 
A = spconvert(A);

% Define the exact solution vector of all ones
n = size(A, 1); % Get the size of the matrix (n x n)
x_exact = ones(n, 1);

%restart parameter fixed to 50
restart = 50;
tol = 1e-8; % Tolerance for GMRES
maxit = 550;

% Define the drop tolerances for the ILU preconditioner
drop_tolerances = [2e-2, 1e-2, 3e-3, 1e-3, 1e-4, 1e-5];

% Compute the right-hand side b
b = A * x_exact;

% Initialize table to store results
results = table('Size', [6 7], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'DropTolerance', 'Iterations', 'tprec', 'tsol', 'TotalTime', 'FinalResidualNorm', 'Density'});

% Initialize figure for convergence profiles
figure;
hold on;

% Loop over each drop tolerance and perform the necessary computations
for i = 1:length(drop_tolerances)
    droptol = drop_tolerances(i);
    
    % Time the computation of the ILU preconditioner
    tic;
    setup.type = 'crout';          % ILU type
    setup.droptol = droptol;       % Drop tolerance for ILU
    [L, U] =  ilu(A, setup);  
    tprec = toc;
    
    % Time the solution process
    tic;
    [x, flag, relres, iter, resvec] = gmres(A, b, restart, tol, maxit, L, U);
    tsol = toc;

    % Compute total number of iterations
    total_iterations = (iter(1) - 1) * restart + iter(2);

    % Density of the preconditioner
    density = (nnz(L) + nnz(U) - n) / nnz(A);
    
    % Store the results in the table
    results(i, :) = {droptol, length(resvec) - 1, tprec, tsol, tprec + tsol, relres, density};
    
    % Plot the convergence profile
    plot(0:total_iterations, log(resvec) , 'DisplayName', ['droptol = ', num2str(droptol)]);
end

% Finalize the convergence plot
xlabel('Iteration number');
ylabel('Relative residual norm');
title('Convergence profile of GMRES with ILU preconditioner');
legend('show');
hold off;

% Display the results table
disp('Results Table:');
disp(results);
