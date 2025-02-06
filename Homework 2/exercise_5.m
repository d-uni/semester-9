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

% right-hand side
b = A * x_exact;

% Parameters for GMRES
tol = 1e-10;       
maxit = 550;      
x0 = zeros(n, 1);  

% Solve using GMRES
[x_gmres, iter, resvec, flag] = mygmres(A, b, tol, maxit, x0);


resvec_full = [norm(b - A * x0); resvec(1:iter)];

semilogy(0:iter, resvec_full, 'r*-', 'LineWidth', 1.5);
xlabel('Iteration Number');
ylabel('Residual Norm (log scale)');
title('GMRES Residual Norm vs Iteration');
grid on;


function [x, iter, resvec, flag] = mygmres(A, b, tol, maxit, x0)

% Initialization
n = length(b);
r0 = b - A * x0;
ro0 = norm(r0, 2);
beta = ro0;
v1 = r0 / beta;
V = zeros(n, maxit + 1);
H = zeros(maxit + 1, maxit);
V(:, 1) = v1;
resvec = zeros(maxit, 1);
flag = 0;

% Check initial residual
if beta < tol * norm(b, 2)
    x = x0;
    iter = 0;
    resvec = beta;
    return;
end

% GMRES Iteration
for k = 1:maxit
    vk_plus_1 = A * V(:, k);

    for j = 1:k
        H(j, k) = V(:, j)' * vk_plus_1;
        vk_plus_1 = vk_plus_1 - H(j, k) * V(:, j);
    end

    H(k + 1, k) = norm(vk_plus_1, 2);
    V(:, k + 1) = vk_plus_1 / H(k + 1, k);
    [Q, R] = qr(H(1:k + 1, 1:k));

    e1 = zeros(k + 1, 1);
    e1(1) = beta;
    g = Q' * e1;
    rok = abs(g(k + 1));

    resvec(k) = rok;

    % Check convergence
    if rok < tol * norm(b, 2)
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
end

