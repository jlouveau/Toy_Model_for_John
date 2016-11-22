%GEVD for mutation
k = -0.001;
sigma = 1.2;
mu = -1.5;

x = linspace(-10,6,1000);
y1 = gevpdf(x, k, sigma, mu);
R = gevrnd(k, sigma, mu, 1, 1000);
figure(); histogram(R, 'Normalization', 'probability', 'BinWidth', 0.5);
hold on; plot(x, y1);
title( ['GEVD with k = ' num2str(k) ', sigma = ' num2str(sigma)  ' and mu = ' num2str(mu)]);

positive = 0;
for n = 1:length(R)
    if R(1,n) > 0
        positive = positive +1;
    end
end

positive = positive/length(R)
