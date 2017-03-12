% SDP-based method for integer convex quadratic minimization
%
% Finds an SDP-based lower bound and an upper bound, along with a feasible 
% point that achieves the upper bound, to the following integer problem:
%
% minimize    f(x) = x^T P x + 2q^T x + r
% subject to  x in Z^n,
%
% where P is positive definite.
%
% Input arguments:
%   P: n-by-n real matrix, positive definite
%   q: n-by-1 real vector
%   r: real number, zero if not provided
%   K: number of samples used in the randomized algorithm, set to K = 3n by
%      default
%
% Outputs:
%   lb: an SDP-based lower bound on the optimal value of the problem
%   ub: an upper bound on the optimal value of the problem
%   xhat: an integer point that achieves f(xhat) = ub
%
% Notes:
%   Requires CVX. If MOSEK solver is available, use the code as it is. If 
%   not, comment out line 51, or change the solver to one that's available
%   to the system.

% Copyright (c) 2015. Jaehyun Park and Stephen Boyd

function [lb, ub, xhat] = ICQM(P, q, r, K)

    % Input sanity check
    assert(nargin >= 2);
    n = size(P, 1);
    if nargin < 3
        r = 0;
    end
    if nargin < 4
        K = 3*n;
    end
    chol(P); % Checks positive definiteness of P
    assert(size(q, 1) == n && size(q, 2) == 1);

    % Change of coordinates
    xcts = -P\q;
    xflr = floor(xcts);
    [q, r] = NEW_ORIGIN(xflr, P, q, r); % Translate by -xflr

    % Use CVX to solve SDP
    M = [P, q; q', 0];
    cvx_solver mosek % Needs CVX professional license
    cvx_begin quiet
        variables X(n, n) x(n)
        expression Z
        Z = [X, x; x', 1];
        minimize sum(sum(M.*Z))
        subject to
            Z == semidefinite(n+1);
            diag(X) >= x;
    cvx_end
    lb = cvx_optval;
    X = (X + X')/2; % Force symmetry

    % Randomized algorithm
    mu = x; Sigma = X - x*x';
    % Using eigenvalue decomposition because chol() complains about
    % near-zero negative eigenvalues
    [V, D] = eig(Sigma);
    A = V*sqrt(max(D, 0));

    ub = 0; xhat = zeros(n, 1);
    for k = 1:K
        x = round(mulrandn_cached(mu, A));
        [x, val] = ONEOPT(x, P, q);
        if ub > val
            ub = val; xhat = x;
        end
    end
    
    % Undo the change of coordinates
    lb = lb + r;
    ub = ub + r;
    xhat = xhat + xflr;
end



% Subroutine for change of coordinates

function [q2, r2] = NEW_ORIGIN(x0, P, q, r)
    q2 = P*x0 + q;
    r2 = x0'*P*x0 + 2*q'*x0 + r;
end


% Subroutine for multivariate normal distribution
%
% Samples from a multivariate normal distribution with mean mu
% and square root matrix of covariance, A (covariance = A*A').
%
% Reference:
% http://en.wikipedia.org/wiki/Multivariate_normal_distribution

function x = mulrandn_cached(mu, A)
    n = size(mu, 1);
    z = randn(n, 1);
    x = mu + A*z;
end


% 1-opt greedy descent subroutine
%
% Starting from x, greedily chooses and optimizes over a single coordinate 
% to reduce the objective x^T P x + 2 q^T x the most. Stops when no single 
% coordinate change can reduce the objective. Returns the 1-opt point and 
% the function value.

function [x, val] = ONEOPT(x, P, q)
    g = 2*(P*x+q);
    v = diag(P);
    iters = 0;
    while true
        iters = iters + 1;
        if v >= abs(g)
            break;
        end
        c = round(-g./(2*v));
        diffs = (c.^2).*v + c.*g;
        [~, i] = min(diffs);
        x(i) = x(i) + c(i);
        g = g + 2*c(i)*P(:, i);
    end
    val = x'*P*x + 2*q'*x;
end
