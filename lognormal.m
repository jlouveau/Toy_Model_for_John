%close all

%lognormal distribution for mutation
offset = 3.0;
sigma = 0.5;
mu = 1.9;
factor = 1;

R = factor * (offset - exp( normrnd(mu, sigma, 1, 100000)));
T = normrnd(mu, sigma, 1, 100000);
figure(); histogram(R, 'Normalization', 'probability', 'Binwidth', 0.2);
title( ['lognormal distribution ' num2str(factor) ' (' num2str(offset) ' - exp( N(' num2str(mu) ', ' num2str(sigma) '))']);

positive = 0;
for n = 1:length(R)
    if R(1,n) > 0
        positive = positive +1;
    end
end

positive = positive/length(R)