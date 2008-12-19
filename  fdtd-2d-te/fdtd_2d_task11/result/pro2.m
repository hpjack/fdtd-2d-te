function [Ex] = pro2;

pars = csvread('pars.txt');
ysize = pars(1);
zsize = pars(2);
ds = pars(7);

rf = pars(10);

Ex = csvread('intenEx.dat');    % ��������� ������������� Ex^2
k = ysize / zsize;

 Z = (1:zsize) .* (ds * 1e6);  % ������ � "������������" ��������
 Y = (1:ysize) .* (ds * 1e6);  % ������ � "������������" ��������


% ������������� n
figure();
eps = csvread('refr_index.dat');
imagesc(Z, -Y, eps);
grid on;
set(gca, 'FontSize', 16)
xlabel('z, mkm');
ylabel('y, mkm');
saveas(gcf, 'epsilon.jpg', 'jpg');

% ���������������, ��� ������ 
ss = get(0,'ScreenSize');
pict_width = floor(1 * ss(3));
pict_height = floor(pict_width * k);
figure('Position', [1,1,pict_width, pict_height]);
caxis([0.1, 1]);   % ���������� �������������
colormap('gray');
surface(Z, Y, Ex); 
imagesc(Z, -Y, Ex);
grid;
set(gca, 'FontSize', 16)
xlabel('z, mkm');
ylabel('y, mkm');
shading flat;
saveas(gcf, 'moving.jpg', 'jpg');


% ������������� ������������� �� ���
figure();
for i = 1:zsize      
    U(i) = Ex(floor(ysize/2), i);
end
plot(Z, U);
grid on;
set(gca, 'FontSize', 16)
xlabel('z, mkm');
ylabel('I, rel.units');
saveas(gcf, 'axis.jpg', 'jpg');


%������������� ������������� � ��������� ���������
figure();

for i = 1:ysize       
    F(i) = Ex(i, floor((9.18e-6 / ds)));
end
plot(Y, F);
grid on;
set(gca, 'FontSize', 16)
xlabel('y, mkm');
ylabel('I, rel.units');
saveas(gcf, 'focus.jpg', 'jpg');


 z_me = floor((9.18e-6)/ ds);
 y1 = floor(3.2e-6 / ds);
 y2 = floor(3.8e-6 / ds);
 
 energy = 0;
 for i=y1:y2
     energy = energy + Ex(i, z_me);
 end
 energy

% y1 = floor(1e-6 / ds);
% y2 = floor(7e-6 / ds);
% 
% z1 = floor(1e-6 / ds);
% z2 = floor(16e-6 / ds);
% 
% for k=z1:z2
%  energy = 0;
%  for i=y1:y2
%      energy = energy + Ex(i, k);
%  end
%  en(k) = energy;
% end
% figure;
% plot(en);