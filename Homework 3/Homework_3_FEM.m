clc; 
clear;

% Uncomment the necessary sections to run simulations for each mesh level.
% The program allows to analyze the solution, generate plots, and compute
% summary convergence tables for different refinement levels.
% To run the summary convergence analysis, all meshes must be processed (~ 7 minuts).


% RUN FOR MESH 0
% Solve the problem using the 'mesh0'. 
% This function computes the solution u, processes data for boundary traces,
% tracking points, and assembles required matrices for analysis.
[coord0, topol0, u0, trace_path0, u_boundary0, time_steps0, u_tracking0, n_tracking_points0, analysis_times0, K_1_bc0, rhs_bc0] = GetSolution('mesh0');

% Uncomment this line to generate plots for the stationary solution and solution along the boundary.
plots_of_the_solution(coord0, u0, analysis_times0, trace_path0, u_boundary0, n_tracking_points0, time_steps0, u_tracking0);

% Uncomment this line to display a table of solution values and the
% associated error w.r.t. Tab. 1
table_of_solution_values(n_tracking_points0, analysis_times0, time_steps0, u_tracking0);


%RUN FOR MESH 1
%[coord1, topol1, u1, trace_path1, u_boundary1, time_steps1, u_tracking1, n_tracking_points1, analysis_times1, K_1_bc1, rhs_bc1] = GetSolution('mesh1');
%plots_of_the_solution(coord1, u1, analysis_times1, trace_path1, u_boundary1, n_tracking_points1, time_steps1, u_tracking1);
%table_of_solution_values(n_tracking_points1, analysis_times1, time_steps1, u_tracking1);

%RUN FOR MESH 2
%[coord2, topol2, u2, trace_path2, u_boundary2, time_steps2, u_tracking2, n_tracking_points2, analysis_times2, K_1_bc2, rhs_bc2] = GetSolution('mesh2');
%plots_of_the_solution(coord2, u2, analysis_times2, trace_path2, u_boundary2, n_tracking_points2, time_steps2, u_tracking2);
%table_of_solution_values(n_tracking_points2, analysis_times2, time_steps2, u_tracking2);

%RUN FOR MESH 3
%[coord3, topol3, u3, trace_path3, u_boundary3, time_steps3, u_tracking3, n_tracking_points3, analysis_times3, K_1_bc3, rhs_bc3] = GetSolution('mesh3');
%plots_of_the_solution(coord3, u3, analysis_times3, trace_path3, u_boundary3, n_tracking_points3, time_steps3, u_tracking3);
%table_of_solution_values(n_tracking_points3, analysis_times3, time_steps3, u_tracking3);

%RUN FOR MESH 4
%[coord4, topol4, u4, trace_path4, u_boundary4, time_steps4, u_tracking4, n_tracking_points4, analysis_times4, K_1_bc4, rhs_bc4] = GetSolution('mesh4');
%plots_of_the_solution(coord4, u4, analysis_times4, trace_path4, u_boundary4, n_tracking_points4, time_steps4, u_tracking4);
%table_of_solution_values(n_tracking_points4, analysis_times4, time_steps4, u_tracking4);


% RUN SUMMARY CONVERGENCE ANALYSIS
% To compute and plot the convergence information, you must run all meshes.

% Semi-logarithmic convergence plots for Cholesky solver
%semi_logarithmic_convergence_plots_Cholesky(K_1_bc0, K_1_bc1, K_1_bc2, K_1_bc3, K_1_bc4, rhs_bc0, rhs_bc1, rhs_bc2, rhs_bc3, rhs_bc4);

% Semi-logarithmic convergence plots for Jacobi solver
%semi_logarithmic_convergence_plots_Jacobi(K_1_bc0, K_1_bc1, K_1_bc2, K_1_bc3, K_1_bc4, rhs_bc0, rhs_bc1, rhs_bc2, rhs_bc3, rhs_bc4);

% Summary table of the FEM convergence analysis
%summary_table_for_the_FEM_convergence(u0, u1, u2, u3, u4,coord0, coord1, coord2, coord3, coord4, topol0, topol1, topol2, topol3, topol4);


% Function to compute and display a summary table for the FEM convergence analysis
function summary_table_for_the_FEM_convergence(u0, u1, u2, u3, u4, coord0, coord1, coord2, coord3, coord4, topol0, topol1, topol2, topol3, topol4)

    % Load reference solution data (x, y coordinates and solution values)
    refData = load('solRef.txt'); 
    xRef = refData(:, 1);       
    yRef = refData(:, 2);       
    uRef = refData(:, 3);  

    coord = {coord0, coord1, coord2, coord3, coord4};                      % Node coordinates for each mesh
    topol = {topol0, topol1, topol2, topol3, topol4};                      % Element connectivity for each mesh
    u = {u0, u1, u2, u3, u4};                                              % Solution vectors for each mesh
    
    % Algorithm 3.3 Interpolation routines (from assignment)
    interpRef = scatteredInterpolant(xRef, yRef, uRef); 
    
    % Initialize arrays
    epsilons = zeros(length(coord), 1);
    h = zeros(length(coord), 1);
    r = zeros(length(coord)-1,1);

    % Loop over each mesh to calculate the error and mesh size
    for i = 1:length(coord)
        % Interpolate the reference solution at the mesh's nodes
        uRefInterp = interpRef(coord{i}(:, 1), coord{i}(:, 2));

        % Compute nodal areas for the mesh
        nodalAreas = computeNodalAreas(coord{i}, topol{i}); 
        error_squared = 0;
        for j = 1:length(coord)
            error_squared = error_squared + nodalAreas(j) / 3 * (u{i}(j) - uRefInterp(j))^2;
        end

        % Compute the root mean squared error
        epsilons(i) = sqrt(error_squared);

        % Calculate the mesh size
        h(i) = calculateMeshSize(coord{i}, topol{i});
    end

    % Compute the convergence ratio 
    for i = 1:length(r)
        r(i) = (h(i+1)/h(i))*(h(i+1)/h(i))*epsilons(i)/epsilons(i+1);
    end

    % Prepare data for the table
    level_labels = arrayfun(@(x) sprintf('Mesh %d', x-1), 1:length(coord), 'UniformOutput', false);
    error_column = array2table(epsilons, 'VariableNames', {'Epsilon'});
    ratio_column = array2table([NaN; r], 'VariableNames', {'Ratio'}); 
    mesh_column = array2table(h, 'VariableNames', {'MeshSize'});

    % Display the convergence summary table
    convergence_table = [table(level_labels(:), 'VariableNames', {'Mesh'}), ...
                     mesh_column, ...
                     error_column, ...
                     ratio_column];

    disp('FEM Convergence Summary Table:');
    disp(convergence_table);
end

% Function to calculate the mesh size
function h = calculateMeshSize(coord, topol)
    n_elements = size(topol, 1);
    max_side_lengths = zeros(n_elements, 1);
    
    for e = 1:n_elements
        nodes = topol(e, :);
        
        x1 = coord(nodes(1), :);
        x2 = coord(nodes(2), :);
        x3 = coord(nodes(3), :);

        % Lengths of the sides of the triangle
        side1 = norm(x1 - x2);
        side2 = norm(x2 - x3);
        side3 = norm(x3 - x1);
        
        % Maximum side length for this element
        max_side_lengths(e) = max([side1, side2, side3]);
    end
   
    h = max(max_side_lengths);
end


% Function to generate semi-logarithmic plots of the convergence of the 
% PCG method with incomplete Cholesky preconditioner for different meshes.
function semi_logarithmic_convergence_plots_Cholesky(K_1_bc0, K_1_bc1, K_1_bc2, K_1_bc3, K_1_bc4, rhs_bc0, rhs_bc1, rhs_bc2, rhs_bc3, rhs_bc4)
    figure;
    hold on;
    colors = ['r', 'g', 'b', 'm', 'c']; 
    K_1_bc = {K_1_bc0, K_1_bc1, K_1_bc2, K_1_bc3, K_1_bc4,};
    rhs_bc = {rhs_bc0, rhs_bc1, rhs_bc2, rhs_bc3, rhs_bc4};
    tol = 1e-8;             
    maxit = 1000; 
    for i = 1:length(K_1_bc)
        L = ichol(K_1_bc{i}); 
        [~, ~, ~, ~, ichol_residuals] = pcg(K_1_bc{i}, rhs_bc{i}, tol, maxit, L, L');
        plot(1:length(ichol_residuals), log(ichol_residuals)/log(ichol_residuals(1)), 'Color', colors(i), 'DisplayName', sprintf('Mesh %d', i-1));
    end

    xlabel('Iteration');
    ylabel('Relative Residual Norm');
    title('PCG Convergence (Incomplete Cholesky preconditioner)');
    grid on;
    legend('show');
    hold off;
end

% Function to generate semi-logarithmic plots of the convergence of the 
% PCG method with an incomplete Jacobi preconditioner for different meshes.
function semi_logarithmic_convergence_plots_Jacobi(K_1_bc0, K_1_bc1, K_1_bc2, K_1_bc3, K_1_bc4, rhs_bc0, rhs_bc1, rhs_bc2, rhs_bc3, rhs_bc4)
    
    figure;
    hold on;
    colors = ['r', 'g', 'b', 'm', 'c']; 
    K_1_bc = {K_1_bc0, K_1_bc1, K_1_bc2, K_1_bc3, K_1_bc4,};
    rhs_bc = {rhs_bc0, rhs_bc1, rhs_bc2, rhs_bc3, rhs_bc4};
    tol = 1e-8;             
    maxit = 1000; 
    for i = 1:length(K_1_bc)
        M_jacobi = diag(diag(K_1_bc{i})); %Jacobi preconditioner
        [~, ~, ~, ~, jacobi_residuals] = pcg(K_1_bc{i}, rhs_bc{i}, tol, maxit, M_jacobi);
        plot(1:length(jacobi_residuals), log(jacobi_residuals)/log(jacobi_residuals(1)), 'Color', colors(i), 'DisplayName', sprintf('Mesh %d', i-1));
    end

    xlabel('Iteration');
    ylabel('Relative Residual Norm');
    title('PCG Convergence (Incomplete Jacobi preconditioner)');
    grid on;
    legend('show');
    hold off;
end

% Function to generate plots of the stationary solution and solution along the boundary. 
function plots_of_the_solution(coord, u, analysis_times, trace_path, u_boundary, n_tracking_points, time_steps, u_tracking)
    
    % Plot 1: Stationary solution (3D scatter plot of u over the domain).
    figure;
    scatter3(coord(:, 1), coord(:, 2), u, 50, u, 'filled'); % Scatter plot for u at nodes.
    colormap(jet);
    xlabel('x asis');
    ylabel('y asis');
    zlabel('u(x, y)');
    title('Stationary solution');
    grid on;
    view(3); 

    % Plot 2: Solution trace along the boundary at different times
    figure;
    hold on;
    for i = 1:length(analysis_times)
        plot(trace_path, u_boundary(:, i), '-o', 'DisplayName', sprintf('t = %.2f', analysis_times(i)));
    end
    xlabel('arc length');
    ylabel('u');
    title('Trace of Solution along boundary at different times');
    grid on;
    legend('show'); 
    hold off;

    % Plot 3: Solution at tracking points over time
    figure;
    hold on;
    for i = 1:n_tracking_points
        plot(time_steps, u_tracking(i, :), 'DisplayName', sprintf('P%d', i));
    end
    xlabel('time');
    ylabel('u');
    title('Solution at tracking points over time');
    grid on;
    legend('show');
    hold off;

end

% Function to create and display a table of solution values at specific tracking points and times
function table_of_solution_values(n_tracking_points, analysis_times, time_steps, u_tracking)
    % Inputs:
    %   - n_tracking_points: Number of tracking points (P1, P2, P3)
    %   - analysis_times: Specific times at which the solution is analyzed (2.5, 5.0, etc.)
    %   - time_steps: Array of all simulation time steps
    %   - u_tracking: Matrix of solution values at the tracking points across all time steps
    %
    % Outputs:
    %   - Displays two tables:
    %       1. Table of solution values at specified times for each tracking point
    %       2. Error table comparing these values with a reference table
 
    % Initialize a matrix to store solution values for the specified analysis times
    u_tracking_table = zeros(n_tracking_points, length(analysis_times));

    % Loop over the specified analysis times
    for j = 1:length(analysis_times)
        % Find the closest time index in the time_steps array
        [~, time_idx] = min(abs(time_steps - analysis_times(j)));

        % Loop over the tracking points and store the solution value at the corresponding time
        for i = 1:n_tracking_points
            u_tracking_table(i, j) = u_tracking(i, time_idx);
        end
    end

    % Define labels for the tracking points (P1, P2, P3) and times (2.5, 5.0, etc.)
    point_labels = {'u(P1)', 'u(P2)', 'u(P3)'};
    time_labels = {'2.5', '5.0', '7.5', '10.0'};
    tracking_table = array2table(u_tracking_table', 'RowNames', time_labels, 'VariableNames', point_labels);

    % Display the solution table
    disp('Solution at different times for points P1, P2 and P3:');
    disp(tracking_table);
    
    % Reference solution values from a predefined table (Tab.1)
    u_tracking_from_Tab1 = [ 0.2434390, 0.0772287, 0.0183037;...
        0.6046775, 0.2751718, 0.0928241;...
        0.7454968, 0.4328630, 0.1716526;...
        0.7751273, 0.4805514, 0.2008722
        ];

    tracking_table = array2table(abs(u_tracking_table' - u_tracking_from_Tab1), 'RowNames', time_labels,'VariableNames', point_labels);

    % Display the error table
    disp('Associated error, measure as difference with respect to values in Tab.1:');
    disp(tracking_table);
end

function [coord, topol, u, trace_path, u_boundary, time_steps, u_tracking, n_tracking_points, analysis_times, K_1_bc, rhs_bc] = GetSolution(meshName)
    % Function to solve the FEM problem for a given mesh
    % Inputs:
    %   - meshName: Name of the mesh folder containing input files
    % Outputs:
    %   - coord: Node coordinates [x, y]
    %   - topol: Element connectivity matrix [node1, node2, node3]
    %   - u: Solution vector at the nodes
    %   - trace_path: Values of the trace data along the boundary
    %   - u_boundary: Solution at boundary nodes at specific times
    %   - time_steps: Time vector for the simulation
    %   - u_tracking: Solution values at tracking points over time
    %   - n_tracking_points: Number of tracking points
    %   - analysis_times: Times at which the solution is analyzed
    %   - K_1_bc, rhs_bc: Stiffness matrix 

    baseDir = fullfile('.', 'Input_files', meshName, meshName);

    % Read the files
    coord = readmatrix(fullfile(baseDir, 'coord.txt'));                    % Node coordinates [x, y]
    topol = readmatrix(fullfile(baseDir, 'topol.txt'));                    % Element connectivity [node1, node2, node3]
    bound = readmatrix(fullfile(baseDir, 'bound.txt'));                    % Boundary conditions [node, value]
    trace_data = load(fullfile(baseDir, 'trace.txt'));                     % Boundary trace
    tracking_points = load(fullfile(baseDir, 'track.txt'));                % Tracking points for solution analysis

    %Parameters
    delta_t = 0.02;                                                        % Time step
    t_max = 10.0;                                                          % Maximum simulation time
    theta = 0.5;                                                           % Theta scheme (0.5 for Crank-Nicholson)
    n_nodes = size(coord, 1);
    n_elements = size(topol, 1);

    % Initialization
    u = zeros(n_nodes, 1);   
    t = 0;         
    tol = 1e-8;             
    maxit = 1000; 

    % Sparse matrices for FEM computations
    H = sparse(n_nodes, n_nodes);                                          % Stiffness matrix
    M = sparse(n_nodes, n_nodes);                                          % Mass matrix
    f = zeros(n_nodes, 1);                                                 % Right-hand side
 
    % Boundary and trace data
    trace_nodes = trace_data(:, 1);                                        % Nodes along the trace
    trace_path = trace_data(:, 2);                                         % Path values corresponding to trace nodes
    analysis_times = [t_max / 4, t_max / 2, 3 * t_max / 4, t_max];         % Analysis times
    u_boundary = zeros(length(trace_nodes), length(analysis_times)); 

    % Tracking points for analysis
    n_tracking_points = length(tracking_points);                           % Number of tracking points
    u_tracking = zeros(n_tracking_points, floor(t_max / delta_t) + 1);     % Solution values at tracking points over time
    time_steps = 0:delta_t:t_max;                                          % Time vector for the simulation

    % Assembly of global matrices (stiffness and mass)
    for e = 1:n_elements
        nodes = topol(e, :);                                               % Nodes of the element (3 шт)
        coords = coord(nodes, :);                                          % Coordinates of these nodes (3х2)
        % Compute element stiffness matrix, mass matrix, and area
        [H_e, M_e, area] = computeElementMatrices(coords);
     
        for i = 1:3
            for j = 1:3
               H(nodes(i), nodes(j)) = H(nodes(i), nodes(j)) + H_e(i, j);
               M(nodes(i), nodes(j)) = M(nodes(i), nodes(j)) + M_e(i, j);
            end
        end
    end

    % Initialize indices for analysis and time-stepping
    current_time_idx = 1;
    time_idx = 1;

    % Time-stepping loop
    while t <= t_max

        % Assemble system matrices for the current time step
        K_1 = M / delta_t + theta * H;
        K_2 = M / delta_t - (1 - theta) * H;
        rhs = K_2 * u;                                                     % Compute the right-hand side vector

        % Apply boundary conditions to the system
        [K_1_bc, rhs_bc] = applyBoundaryConditions(K_1, rhs, bound, t, t_max);

        % Solve the system using PCG method
        L = ichol(K_1_bc);
        [u, flag, relres, iter, resvec] = pcg(K_1_bc, rhs_bc, tol, maxit, L, L');
    
        t = t + delta_t;
    
        %For Analysis 1 (store the solution at boundary nodes)
        if current_time_idx <= length(analysis_times) && abs(t - analysis_times(current_time_idx)) < delta_t / 2
            u_boundary(:, current_time_idx) = u(trace_nodes);
            current_time_idx = current_time_idx + 1;
        end
    
        %For Analysis 2 (store the solution at tracking points)
        for i = 1:n_tracking_points
            u_tracking(i, time_idx) = u(tracking_points(i));
        end
        time_idx = time_idx + 1;
    
    end

end



% Function to compute nodal areas for all nodes in the mesh
function nodalAreas = computeNodalAreas(coord, topol)
    % Inputs:
    %   - coord: Matrix of node coordinates [x, y]
    %   - topol: Element connectivity matrix [node1, node2, node3]
    % Outputs:
    %   - nodalAreas: Array containing the area contribution of each element to the nodes

    n_nodes = size(coord, 1);
    nodalAreas = zeros(n_nodes, 1);
    
    for e = 1:size(topol, 1)                                               % Loop over each element in the mesh
        nodes = topol(e, :);                                               % Nodes belonging to the current element
        coords = coord(nodes, :);                                          % Coordinates of these nodes
        x = coords(:, 1);                                                  % x-coordinates of the nodes
        y = coords(:, 2);                                                  % y-coordinates of the nodes

        % Compute the area of the triangular element using the determinant formula
        area = 0.5*abs(det([1, x(1), y(1); 1, x(2), y(2); 1, x(3), y(3)]));
        nodalAreas(nodes) = nodalAreas(nodes) + area;
    end
end

% Function to compute the stiffness and mass matrices for a single triangular element
function [H_e, M_e, area] = computeElementMatrices(coords)
    % Inputs:
    %   - coords: Matrix of node coordinates for the element [x1, y1; x2, y2; x3, y3]
    % Outputs:
    %   - H_e: Element stiffness matrix (3x3)
    %   - M_e: Element mass matrix (3x3)
    %   - area: Area of the triangular element

    x = coords(:, 1); 
    y = coords(:, 2);
    area =  0.5*abs(det([1, x(1), y(1); 1, x(2), y(2); 1, x(3), y(3)]));

    b = [y(2) - y(3); y(3) - y(1); y(1) - y(2)];
    c = [x(3) - x(2); x(1) - x(3); x(2) - x(1)];

    H_e = (b * b' + c * c') /(4*area);
    M_e = area / 12 * [2, 1, 1; 1, 2, 1; 1, 1, 2];
end

% Function to apply boundary conditions to the global stiffness matrix and rhs vector
function [K, rhs] = applyBoundaryConditions(K, rhs, bound, t, t_max)
    % Inputs:
    %   - K: Global stiffness matrix
    %   - rhs: Global right-hand side vector
    %   - bound: Boundary conditions matrix [node, value]
    %   - t: Current simulation time
    %   - t_max: Maximum simulation time
    % Outputs:
    %   - K: Modified stiffness matrix with boundary conditions applied
    %   - rhs: Modified RHS vector with boundary conditions applied

    for i = 1:size(bound, 1)                                               % Loop over all boundary nodes
        node = bound(i, 1);                                                % Node index
        value  = bound(i, 2);                                              % Prescribed boundary value

        % Modify the boundary value 
        if t <= t_max / 2
          value = (2 * t / t_max) * value;
        else
          value = value;
        end

        % Modify the matrix and RHS
        % Algorithm 3.2 (from assignment) Enforcement of Dirichlet boundary condition
        K(node, :) = 0;                                                   
        rhs = rhs - K(:, node) * value;                                   
        rhs(node) = value;
        K(:, node) = 0;
        K(node, node) = 1;
    end
end
