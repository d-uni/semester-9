clc;
clear;

A = delsq(numgrid('S', 102)); 
L = ichol(A); % Cholesky preconditioner
n = size(A, 1);
b = A * ones(n, 1); 
tol = 1e-8; 
maxit = 50;

% Run mypcg
[x_mypcg, resvec_mypcg, iter_mypcg] = mypcg(A, b, tol, maxit, L);

% Run MATLAB's PCG for comparison
[x_pcg, flag, relres, iter_pcg, resvec_pcg] = pcg(A, b, tol, maxit, L, L');

fprintf('Iterations (mypcg): %d\n', iter_mypcg);
fprintf('Iterations (pcg): %d\n', iter_pcg);

pcg_iters = length(resvec_pcg) - 1;
mypcg_iters = length(resvec_mypcg) - 1;

semilogy(0:pcg_iters, resvec_pcg, 'r*-', 0:mypcg_iters, resvec_mypcg, 'go-');
legend('MATLAB pcg', 'mypcg');
xlabel('Iteration number');
ylabel('Residual norm');
title('Comparison of PCG Convergence');


function [x, resvec, iter] = mypcg(A, b, tol, maxit, L)

    n = length(b);
    x = zeros(n, 1);                
    r = b - A * x;                   
    M_inv_r = L \ (L' \ r);          
    p = M_inv_r;                     
    rho_old = r' * M_inv_r;          

    resvec = zeros(maxit + 1, 1);    
    resvec(1) = norm(r);             
    iter = 0;

    % PCG loop
    while norm(r) > tol * norm(b) && iter < maxit
        iter = iter + 1;
        z = A * p;
        alpha = rho_old / (p' * z);
        x = x + alpha * p;
        r = r - alpha * z;
        M_inv_r_new = L \ (L' \ r);
        rho_new = r' * M_inv_r_new;
        beta = rho_new / rho_old;
        p = M_inv_r_new + beta * p;

        rho_old = rho_new;
        resvec(iter + 1) = norm(r);  
    end
    resvec = resvec(1:iter + 1);
end