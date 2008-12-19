function [Ex] = pro2;

pars = csvread('pars.txt');
ds = pars(7);

d = 0.25e-6;
L = 3e-6;
n = 2.83;
N = floor(2.5e-6 / ds);
for i=1:N
    t1 = d * n  / (2 * (n - 1));
    t2 = cosh(pi * i * ds / 2 / L);
    r(i)= t1 * (1 - 1 / t2);
end
x = (1:N).*(ds*1e6);
R = r.*1e6;

figure;
set(gca, 'FontSize', 16)
plot(x, R);
xlabel('y, mkm');
ylabel('r, mkm');

M = floor(2.5e-6 / d);
for i=1:M
    ind2 = floor((i)*d * N / 2.5e-6);
    
    re(i) = floor(r(1 + ind2) / ds + 0.5) * ds;
end

re = re.*1e6;

figure;
set(gca, 'FontSize', 16)
x2 = (1:M) * d * 1e6;
plot(x2, re, x, R);
xlabel('y, mkm');
ylabel('r, mkm');