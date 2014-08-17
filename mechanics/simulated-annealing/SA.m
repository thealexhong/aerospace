% This MATLAB code is the simulated annealing (SA) algorithm

% SA.m
% This function is a probabilistic search algorithm analogous to annealing
% of metals. It can solve bound constrained optimization problems.

function [xopt, fopt, bFeasible] = SA (x0, lb, ub, epsilon, ...
                            maxiter, Tstart, c, opt, penalty, P, bGraph, fig)
% Inputs
% x0: inital guess solution
% lb: vector containing lower bounds on the design variables
% ub: vector containing upper bounds on the design variables
% epsilon: step-size controlling magnitude of perturbation to design
... variables
% maxiter: maximum number of iterations
% Tstart: starting temperature
% c: cooling schedule parameter
% opt: passed onto objfcn for evaluating chosen objective function
% penalty: passed onto objfcn for choosing penalty type
% P: penalty parameter
% bGraph: a boolean to display convergence trend
% fig: which figure should the program display it on
%
% Outputs
% xopt: the optimal solution
% fopt: the optimal function output given the optimal solution

% Checks if initial guess lies within bound
if (~isempty(find(x0 < lb, 1)) || ~isempty(find(x0 > ub, 1)))
    error ('The initial guess lies outside of the defined bounds');
end

x = x0; % set initial guess
T = Tstart; % initialize temperature
t = 0; % initial time
[fopt, ~] = objfcn(x, opt, penalty, P); % best is initialize to start
xopt = x;
trend = zeros(1, maxiter); % for graphing

while t < maxiter
    T = schedule(T, c, t); % Update temperature every iteration    
    % Check if temperature has reached termination
    if (T == 0)
        break;
    end
    
    x_ = move (x, lb, ub, epsilon);
    [f_, ~] = objfcn(x, opt, penalty, P);
    [f__, ~] = objfcn(x_, opt, penalty, P);
    deltaE =  f_ - f__;
    if (deltaE > 0)
        x = x_;
    else
        if (exp(deltaE / T) >= rand())
            x = x_;
        end
    end
    
    % Best Solution: Backed up every iteration
    [f, bFeasible] = objfcn(x, opt, penalty, P);
    if (f < fopt);
        xopt = x;
        fopt = f;
    end
    trend(t+1) = f;
    t = t + 1;
end

% Graph the trend
if (bGraph)
    trend = trend(1:t); % useful part
    figure(fig);
    plot (1:t, trend);
    title(sprintf('Convergence trend for %d iterations completed', t));
    xlabel('Iterations');
    ylabel('Objective function');
end