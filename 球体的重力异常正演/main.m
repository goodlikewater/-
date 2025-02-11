clear all;
G=6.67*1e-2;
R=50;
D=100;
sigma=1;
x=-200:10:200;
y=-200:10:200;
M=(4/3)*pi*(R-3)*sigma;
for i=1:size(x,2)
    for j=1:size(x,2)
        g(j,i)=(G*M*D)/((x(i)^2+y(j)^2+D^2)^1.5);
    end
end
%等值线图
subplot(121);
contourf(x,y,g,12);
colorbar;
box on;
xlabel('X/m');
ylabel('Y/m');
%三维曲面图
subplot(122);
surf(x,y,g);
colorbar;
xlabel('X/m');
ylabel('Y/m');
zlabel('Deltag/g.u.');
