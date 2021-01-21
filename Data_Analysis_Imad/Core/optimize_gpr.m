function [hyperpar,lnL] = optimize_gpr( X, Y, theta0 )

n = size(X,1);
problem.objective = @(theta) gpr_likelihood( X, Y, theta );
problem.solver = 'fmincon';
problem.lb = zeros(n+2,1);

if nargin == 3
    problem.x0 = theta0;
else
    problem.x0 = 1e-3*ones(1,n+2);
end

problem.options = optimoptions('fmincon');
problem.options.Display = 'iter';
% problem.options.PlotFcns = @optimplotfval;
% problem.options.CheckGradients = true;
problem.options.UseParallel = true;
problem.options.SpecifyObjectiveGradient = true;

[hyperpar,lnL] = fmincon(problem);