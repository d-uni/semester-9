clc;
clear;

A = gallery('wathen', 100, 100);
n = size(A, 1); % Size of the matrix

x_exact = rand(n, 1);
b = A * x_exact;

% Parameters for PCG
tol = 1e-8; 
maxit = 1000; 

[x_noprec,~, resvec_noprec, iter_noprec, ~] = pcg(A, b, tol, maxit);
M_jacobi = diag(diag(A));
[x_jacobi,~, resvec_jacobi, iter_jacobi, ~] = pcg(A, b, tol, maxit, M_jacobi);

% Solve with IC(0) preconditioner
L_ic0 = ichol(A); 
[x_ic0, ~, resvec_ic0, iter_ic0, ~] = pcg(A, b, tol, maxit, L_ic0, L_ic0');

fprintf('PCG Iterations:\n');
fprintf('Non-preconditioned: %d\n', iter_noprec);
fprintf('Jacobi:             %d\n', iter_jacobi);
fprintf('IC(0):              %d\n', iter_ic0);


