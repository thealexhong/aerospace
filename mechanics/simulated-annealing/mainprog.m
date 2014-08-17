% This MATLAB code implements and analyzes the simulated annealing
... algorithm applied on a bump test function and 10 truss structure

% mainprog.m
% This code is the main function that controls the flow of the program


function mainprog (opt)

% clear
clc;
close all;
fclose all;

switch opt
    case 1
        
        % Get Tstart, lb, ub, maxiter, c, and epsilon
        bumpInput();
        
        % Set testing parameters
        nRuns = 50; % number of runs to determine convergence trend
        c = 0.8:0.01:0.99; % the cooling schedule parameter to test
        % overwrite epsilon defined in bumpInput with test values
        epsilon = [0.3 0.5 0.7 0.9 0.99];
        
        for i = 1:size(epsilon, 2) % find the best epsilon
            fopt_avg = zeros(1, size(c, 2));
            for j = 1:size(c, 2) % find the best c
                for k = 1:nRuns
                    % Random initial guess within range
                    x0 = (ub - lb) .* rand() + lb;
                    [~, fopt, ~] = ...
                        SA(x0, lb, ub, epsilon(i), ...
                            maxiter, Tstart, c(j), opt, 'onepass', ...
                            1000000, 0, 2); % if graphed, it's on figure 2
                    fopt_avg(j) = ...
                        fopt_avg(j) + fopt; % Divide by total after
                end
            end
        
            % Display convergence trend of the test runs
            figure (1)
            plot (c, fopt_avg / nRuns); % average optimal function value
            hold all;
        end
        
        hold off;
        xlabel('Cooling parameter value, c');
        ylabel('Average objective function value');
        title(sprintf(['Average optimum function values' ...
            ' (over %d runs each) vs. cooling parameter'], nRuns));
        legend (num2str(epsilon'));
        
    case 2 
        
        trussInput();
        
        % Set testing parameters
        % NOTE: Increasing this value will take a lot of computational
        % time.
        nRuns = 5;
        
        % Perform SA multiple times to find the average minimum weight of
        % the structure
        % Measure the performace of the one-pass penalty method
        fopts = zeros(1, nRuns);
        for i = 1:nRuns    
            x0 = (ub - lb) .* rand() + lb;
            [~, fopt, bFeasible] = ...
                        SA(x0, lb, ub, epsilon,...
                        maxiter, Tstart, c, opt, 'onepass', 1000000, 1, 1);
            fopts(i) = fopt;
        end
        fprintf (['The avarage solution with one-pass penalty '...
            'is %f kg with\nstd dev %f kg/10 over %d runs.\n\n'],...
            mean(fopts), std(fopts), nRuns);
        fprintf (['The minimum solution with one-pass penalty '...
            'is %f kg/10 over %d runs.\n'],...
            min(fopts), nRuns);
        if (bFeasible)
            disp('== Solution satisfies all constraints (Feasible!) ==');
        end

        
        % Measure impact of penalty parameter P on the convergence trends
        % of quadratic penalty function approach
        P = [0.1, 100, 1000];
        fopts_ = zeros(size(P, 2), nRuns);
        for i = 1:size(P, 2)
            for j = 1:nRuns
                x0 = (ub - lb) .* rand() + lb;
                [~, fopt, ~] = SA(x0, lb, ub, epsilon,...
                        maxiter, Tstart, c, opt, 'quadratic', P(i), 1, i + 1);
                fopts_(i, j) = fopt;
            end
        end
        fprintf(['Quadratic Penalty Approach Statistics over %d runs\n\n' ...
                 '      Penalty    |    Mean   |    StdDev\n'], nRuns);
             format shortG;
        disp([P', mean(fopts_,2), std(fopts_, 0, 2)]);
        
        
    otherwise
        error (sprintf(...
               ['Invalid input for maindriver.\n'...
                'Please use the following inputs for maindriver(opt):\n'...
                'opt = 1: SA on bump function for finding best c, epsilon\n'...
                'opt = 2: SA on 10 beam truss structure\n']));
end

