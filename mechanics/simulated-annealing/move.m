% Normalizes the design variable
... vector, perturbs it, and then de-normalizes it.
% Added Features:
% Normalize design variable vector x before perturbation
% De-normalize design variable vector x after perturbation

function [z] = move(x, lb, ub, epsilon)
%
% [z] = move(x, lb, ub, epsilon)
%
% This matlab function randomly perturbs the design variable
% x such that the perturbed vector satisfies the bound 
% constraints. 
%
% Inputs: 
% ------
% x :      Design variable vector
% lb:      Vector containing lower bounds on the design variables
% ub:      Vector containing upper bounds on the design variables
% epsilon: A parameter controlling the magnitude of perturbation.
%          It is recommended that this parameter is set to a value 
%          between 0.1 and 0.3, if the design variables are normalized to
%          [0,1]. Try some typical values and see what impact 
%          this parameter makes on the convergence trends.
%
% Output:
% ------
% z:  perturbed design variable vector satisfying the bound constraints
%
%

  x = (x - lb) ./ (ub - lb); % Normalize design variable
  
  n = length(x);  % Extract the number of design variables
  flag = 0;
  while flag == 0
        ind = ceil(rand*n); % randomly generate an integer between 1 and n
                            % to select the design variable that will be 
                            % perturbed to generate the move
        z = x;
        z(ind) = x(ind) + epsilon*(-1 + rand*2);   % Perturb the random chosenly 
                                                   % design variable by epsilon*U[-1,1]
        if (z(ind) < 0 || z(ind) > 1)              % If bound constraints are violated
           flag = 0;                               % generate a new perturbation
        else                                       % Good job, let's return
           flag = 1;
        end
  end
  
  z = (ub - lb) .* z + lb; % De-normalize design variable

