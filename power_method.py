%% ========== POWER METHOD WITH ITERATION DISPLAY ==========
% This script implements the Power Method and displays all iterations

% Define the matrix
A = [2   2   2;
     2/3   5/3   5/3;
     1   5/2  11/2];

% Algorithm parameters
epsilon = 1e-9;
maxIterations = 10;

%% First initial guess: x0 = [1; 0; 0]
fprintf('\n=== First Initial Guess: [1; 0; 0] ===\n');
x0 = [1; 0; 0];

% Initialize variables
x = x0 / norm(x0, inf);
iter = 0;
converged = false;

while ~converged && iter < maxIterations
    iter = iter + 1;

    % Power Method steps
    y = A * x;
    [max_val, idx] = max(abs(y));
    mu = y(idx);
    x_new = y / mu;

    % Display current iteration
    fprintf('\nIteration %d:\n', iter);
    fprintf('mu = %.6f\n', mu);
    fprintf('x = [%.6f, %.6f, %.6f]\n', x_new(1), x_new(2), x_new(3));

    % Check convergence
    if norm(x_new - x, inf) < epsilon
        converged = true;
    end

    x = x_new;
end

fprintf('\nConverged in %d iterations!\n', iter);
fprintf('Dominant eigenvalue: %.6f\n', mu);
fprintf('Eigenvector: [%.6f, %.6f, %.6f]\n', x(1), x(2), x(3));
