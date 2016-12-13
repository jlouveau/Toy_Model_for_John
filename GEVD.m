close all
%GEVD for mutation
k = -0.6;
sigma = 1.8;
mu = -1.9;

x = linspace(-10,3,1000);
y1 = gevpdf(x, k, sigma, mu);
R = gevrnd(k, sigma, mu, 1, 100000);
figure(); histogram(R, 'Normalization', 'probability', 'Binwidth', 0.2);
%hold on; plot(x, y1);
%title( ['GEVD with k = ' num2str(k) ', sigma = ' num2str(sigma)  ' and mu = ' num2str(mu)]);

positive = 0;
for n = 1:length(R)
    if R(1,n) > 0
        positive = positive +1;
    end
end

positive = positive/length(R)
