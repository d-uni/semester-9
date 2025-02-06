clc;
clear;

data = load('mat13041.rig'); 

rows = data(:, 1);       
cols = data(:, 2);       
vals = data(:, 3);       

n = max(max(rows), max(cols));
A = sparse(rows, cols, vals, n, n);

n = size(A, 1);


x_exact = 1 ./ sqrt((1:n)'); 
b = A * x_exact;

% Compute the ILU(0.1) preconditioner
setup.type = 'crout';        
setup.droptol = 0.1;        
[L, U] = ilu(A, setup);

% Parameters for GMRES
tol = 1e-10;       
maxit = 550;       
x0 = zeros(n, 1); 

% My Preconditioned GMRES Implementation
[x_precgmres, iter_prec, resvec_prec, flag_prec] = myprecgmres(A, b, tol, maxit, x0, L, U);
true_residual_prec = norm(b - A * x_precgmres);

% Matlab Built-in GMRES with Preconditioning
[x_gmres_builtin, flag_builtin, relres_builtin, iter_builtin, resvec_builtin] = gmres(A, b, [], tol, maxit, L, U, x0);
true_residual_builtin = norm(b - A * x_gmres_builtin);


results_table = table(...
    [iter_prec; iter_builtin(2)], ...                  
    [flag_prec; flag_builtin], ...                   
    [resvec_prec(end); resvec_builtin(end)], ...      
    [true_residual_prec; true_residual_builtin], ...  
    'VariableNames', {'Iterations', 'Flag', 'PreconditionedResidual', 'TrueResidualNorm'}, ...
    'RowNames', {'MyGMRES', 'MatlabGMRES'});


disp(results_table);




function [x, iter, resvec, flag] = myprecgmres(A, b, tol, maxit, x0, L, U)

% Initialize
n = length(b);
r0 = b - A * x0;
z0 = U \ (L \ r0);  % Solve (LU)^-1 * r0
beta = norm(z0, 2);
V = zeros(n, maxit + 1);
H = zeros(maxit + 1, maxit);
V(:, 1) = z0 / beta;
resvec = zeros(maxit, 1);
flag = 0;

% Preconditioned norm of b
z_b = U \ (L \ b);
norm_z_b = norm(z_b, 2);

% Check initial residual
if beta < tol * norm_z_b
    x = x0;
    iter = 0;
    resvec(1) = beta;
    return;
end

% GMRES Iteration
for k = 1:maxit
   
    vk_plus_1 = A * V(:, k);
    z = U \ (L \ vk_plus_1);  

    for j = 1:k
        H(j, k) = V(:, j)' * z;
        z = z - H(j, k) * V(:, j);
    end

    H(k + 1, k) = norm(z, 2);
    V(:, k + 1) = z / H(k + 1, k);
    [Q, R] = qr(H(1:k + 1, 1:k));
    e1 = zeros(k + 1, 1);
    e1(1) = beta;
    g = Q' * e1;
    rok = abs(g(k + 1));

    resvec(k) = rok;

    % Check convergence
    if rok < tol * norm_z_b
        y = R(1:k, 1:k) \ g(1:k);
        x = x0 + V(:, 1:k) * y;
        iter = k;
        resvec = resvec(1:k);
        return;
    end
end

y = R(1:maxit, 1:maxit) \ g(1:maxit);
x = x0 + V(:, 1:maxit) * y;
iter = maxit;
resvec = resvec(1:maxit);
flag = 1; 
end
