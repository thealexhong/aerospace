% objfcn.m
% This function returns the objective function given the design variable
... vector x. A simple one-pass penalty function is added to the
... objective function value if any of the constraints are violated.

% The objective function is programmed to be the bump test function in
... the assignment outline.

function [f, bFeasible] = objfcn(x, opt, penalty, P)
% Variables
% x: design variable vector x
% f: objective function
% opt: depending on the question, which objective function to use
% penalty: specifies which penalty method to use
% P: penalty parameter
% f: evaluated objective function
% bFeasible: boolean to specified whether the solution is within constraint
... and bound

bFeasible = true;

switch opt
    case 1
        n = 2;
        % Bump test function
        f = -abs( sum(cos(x) .^ 4) - 2 * prod(cos(x) .^ 2)) /...
            sqrt( sum((1:n) * (x .^ 2)') );

        % A simple one-pass penalty function
        % Add a large number to the objective function value if any of the
        ... constraints are violated.
        if (prod(x) <= 0.75 || sum(x) >= (15 * n / 2))
            f = f + P;
            bFeasible = false;
        end
    case 2
        
        % Aluminum Properties
        % density of aluminium in kg/m^3 (2.70 g/cm^3)
        rho = 2.70 * (1 / 1000) * (100 / 1)^3;
        yieldstr_tensile = 276E6; % maximum tensile stress (Pa)
        yieldstr_compr = -yieldstr_tensile; % maximum compressive stress (Pa)
        
        area = x; % design variable
        % Store element data from Assignment 1
        addpath Assignment1_Pkg/
        [elementData] = maindriver(2, area); % Use assignment 1
        rmpath Assignment1_Pkg/
        
        f = 0;
        % calculate total volume (m^3)
        for i = 1:size(elementData, 2)
            f = f + elementData(i).A * elementData(i).L;
        end
        f = f * rho; % calculate total mass (kg) (constant density)
        
        switch penalty
            case 'onepass'
                % one-pass penalty function method
                for i = 1:size(elementData, 2)
                    sigma = elementData(i).stress;
                    if (sigma > yieldstr_tensile || sigma < yieldstr_compr)
                        f = f + P;
                        bFeasible = false;
                        break;
                    end
                end
            case 'quadratic'
                % quadratic penalty function method
                q = 0;
                for i = 1:size(elementData, 2)
                    sigma = elementData(i).stress;
                    q = q + (max(sigma - yieldstr_tensile, 0))^2 +...
                        (max(yieldstr_compr - sigma, 0))^2;
                end
                if (q > 0)
                    bFeasible = false;
                end
                f = f + P * q;
            otherwise
                error (['Invalid parameter for penalty.'...
                    ' Use ''onepass'' or ''quadratic''.']);
        end
    otherwise
        error ('Invalid parameter for opt in objfcn (x, opt).');
end