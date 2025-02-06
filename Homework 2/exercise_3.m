clc;
clear;
function [x, resvec, iter] = mypcg(A, b, tol, maxit, L)

    % Initialization
    n = length(b);
    x = zeros(n, 1);                 
    r = b - A * x;                   
    M_inv_r = L \ (L' \ r);          
    p = M_inv_r;                     
    rho_old = r' * M_inv_r;         

    resvec = zeros(maxit + 1, 1);    
    resvec(1) = norm(r);             
    iter = 0;

    %  PCG loop
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
% Parameters
n = 1e4; % Matrix size
tol = 1e-8; % Convergence tolerance
maxit = 50; % Maximum number of iterations

A1 = diag([200*(1:5), ones(1, n-5)]);
b = rand(n, 1);
M = eye(n);
[x_mypcg, resvec, iter] = mypcg(A1, b, tol, maxit, M);
fprintf('PCG completed in %d iterations.\n', iter);

semilogy(0:iter, resvec, 'o-', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('Residual Norm (log scale)');
title('PCG Residual Norm vs Iteration Number');
grid on;





