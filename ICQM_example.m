% Example script for running ICQM

% Construct integer least squares problem
% minimize    ||Ax - b||^2
% subject to  x in Z^n

n = 100;
A = randn(2*n, n) / 10; b = randn(2*n, 1);

% ||Ax - b||^2 = x^T P x + 2q^T x + r
P = A'*A;
q = -A'*b;
r = b'*b;

tic;
[lb, ub, xhat] = ICQM(P, q, r);
elapsed = toc;

fprintf('Elapsed time: %.5f seconds\n', elapsed);
fprintf('Optimal value is between %.5f and %.5f\n', lb, ub);
fprintf('2-norm of the suboptimal solution: %.5f\n', norm(xhat));
