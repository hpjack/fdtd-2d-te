function [ output_args ] = symmetry( input_args )

eps = csvread('refr_index.dat');
ep1 = eps(1:226, :);
ep2 = eps(227:452, :);
n = size(ep2)
for i=1:n(1)
    ep3(i, :) = ep2(n(1)+1-i,:);
end

figure
imagesc(((ep1)./8.85418782e-12))
figure
imagesc(((ep1 - ep3)./8.85418782e-12))