clc;
clear;

% Parameters
nx_values = [102, 202, 402, 802]; 
maxit = 5000; 
tol = 1e-8; 

results = zeros(length(nx_values), 5); % Columns: h, CG, PCG (IC0, IC10^-2, IC10^-3)

for idx = 1:length(nx_values)
    nx = nx_values(idx);
    A = delsq(numgrid('S', nx));
    n = size(A, 1);

    h = 1 / (nx - 1);

    % Right-hand side b for exact solution x
    x_exact = 1 ./ sqrt((1:n)');
    b = A * x_exact;

    % Solve with non-preconditioned CG
    [~, ~, ~, iter_cg, ~] = pcg(A, b, tol, maxit);

    % Preconditioners
    L_IC0 = ichol(A); 
    opts1.type = 'ict'; opts1.droptol = 1e-2;
    L_IC2 = ichol(A, opts1); 
    opts2.type = 'ict'; opts2.droptol = 1e-3;
    L_IC3 = ichol(A, opts2); 

    [~, ~, ~, iter_pcg_ic0, ~] = pcg(A, b, tol, maxit, L_IC0, L_IC0');

    [~, ~, ~, iter_pcg_ic2, ~] = pcg(A, b, tol, maxit, L_IC2, L_IC2');

    [~, ~, ~, iter_pcg_ic3, ~] = pcg(A, b, tol, maxit, L_IC3, L_IC3');

    results(idx, :) = [h, iter_cg, iter_pcg_ic0, iter_pcg_ic2, iter_pcg_ic3];
end

fprintf('  h         CG   PCG_IC0   PCG_IC10^-2   PCG_IC10^-3\n');
for idx = 1:length(nx_values)
    fprintf('%8.5f   %4d   %9d   %12d   %12d\n', results(idx, :));
end

fprintf('\n \n  1/h  CG_Iterations/h^-1  PCG_IC0_Iterations/h^-1  PCG_IC10^-2_Iterations/h^-1   PCG_IC10^-3_Iterations/h^-1\n');
for idx = 1:length(nx_values)
    fprintf('%8.5f   %4d   %9d   %12d   %12d\n', 1/results(idx, 1)...
        , results(idx, 2)*results(idx, 1)...
        , results(idx, 3)*results(idx, 1) ...
        , results(idx, 4)*results(idx, 1)...
        , results(idx, 5)*results(idx, 1));
end
