function [smpl,lnP,lnL,accept] = optimize_gpr_bayes( X, Y, theta0, nsamples )

theta0 = theta0(:)';                                                        % Initial value

lnpdf = @(theta)gpr_posterior(X,Y,theta);                                   % Posterior distribution
proprnd = @(x)proprand(x);                                                  % Proposal random number generator

if nargin == 3                                                              % Number of samples
    nsamples = 1e3;
end

% Metropolis-Hastings algorithm to produce samples from the posterior
fprintf('Sampling the GPR posterior with the Metropolis-Hasting algorithm\n')
[smpl,lnP,lnL,accept] = my_mhsample(theta0,nsamples,lnpdf,proprnd);
fprintf('Done!!!\n')

%-- Posterior probability -------------------------------------------------
function [lnP,lnL] = gpr_posterior(X,Y,theta)

% Fetching information from input
theta = theta(:);                                                           % Parameters in linear scale
n = size(X,1);                                                              % Number of parameters
L = log(theta(3:end));                                                      % Log-scale parameters
ln_sigmaU2 = log(theta(1));                                                 % Log-noise variance
ln_sigmaF2 = log(theta(2));                                                 % Kernel variance

% Parameters of the prior of the scale parameters
muL = zeros(n,1);                                                           % Mean of the log-scale parameters
SigmaL = 1e2*eye(n);                                                        % Covariance of the log-scale parameters

% Parameters of the prior of the variance parameter
muS = -1;                                                                   % Mean of the log-variance
sigmaS = 1e2;                                                               % Variance of the log-variance

% Parameters of the prior of the variance parameter
muS1 = -1;                                                                  % Mean of the log-variance
sigmaS1 = 1e2;                                                              % Variance of the log-variance

% Calculate prior log-probabilities
ln_priorL = log( mvnpdf( L, muL, SigmaL ) );
ln_priorSigma = log( normpdf( ln_sigmaU2, muS, sigmaS ) );
ln_priorSigma1 = log( normpdf( ln_sigmaF2, muS1, sigmaS1 ) );

% Calculate log-likelihood
lnL = gpr_likelihood( X, Y, theta );

% Calculate posterior
lnP = -lnL + ln_priorL + ln_priorSigma + ln_priorSigma1;

%-- Proposal random number generator --------------------------------------
function xnew = proprand( x )

% Fetching information from input
n = size(x,2);                                                              % Size of the parameter vector
Sigma = 1e-1*eye(n);                                                        % Covariance matrix

xnew = exp( log(x) + randn(1,n)*Sigma );                                    % Create a new sample point

%-- Metropolis-Hastings sampling ------------------------------------------
function [x,f,g,accept] = my_mhsample( initial, nsamples, lnpdf, proprnd )

n = length(initial(:));
x = zeros(nsamples,n);
f = zeros(nsamples,1);
g = zeros(nsamples,1);
accept = zeros(nsamples,1);

x(1,:) = initial(:)';
[f(1),g(1)] = lnpdf(x(1,:));

for i=2:nsamples
    
    % Draw a new sample
    x_new = proprnd(x(i-1,:));
    
    % Calculate the acceptance ratio
    [f_new,g_new] = lnpdf(x_new);
    accept(i) = min(0,f_new-f(i-1));
    
    % Accept / reject
    u = log(rand);
    if accept(i) >= u
        x(i,:) = x_new;
        f(i) = f_new;
        g(i) = g_new;
    else
        x(i,:) = x(i-1,:);
        f(i) = f(i-1);
        g(i) = g(i-1);
    end
    
end

