x=-40:1:40;
y=-40:1:40;
R=15; %������������(m)
B=50000; %���شų�(nT)
k=0.015; %�Ż��ʣ�SID
r=10; %����뾶(m)
vol=4*pi*r*r*r/3; %�������
mu0=4*pi*1e-7; %������
m=k*B/mu0; %����Ż�ǿ��
M=m*vol; %����ž�
A=45*pi/180;
I=0;

for i=1:size(x,2)
    for j=1:size(y,2)
        Hax(j,i)=(mu0/(4*pi))*M*((2*x(i)*x(i)-y(j)*y(j)-R*R)*cos(I)*cos(A)...
            -3*R*x(i)*sin(I)+3*x(i)*y(j)*cos(I)*sin(A))/((R^2+x(i)^2+y(j)^2)^2.5);
        Hay(j,i)=(mu0/(4*pi))* M *((2*y(j)*y(j)-x(i)*x(i)-R*R)*cos(I)*cos(A)...
            -3*R*y(j)*sin(I)+3*x(i)*y(j)*cos(I)*cos(A))/((R^2+x(i)^2+y(j)^2)^2.5);
        Za(j,i)=(mu0/(4*pi))* M *((2*R*R-x(i)*x(i)-y(j)*y(j))*sin(I)...
            -3*R*x(i)*cos(I)*cos(A)-3*R*y(j)*cos(I)*sin(A))/((R^2+x(i)^2+y(j)^2)^2.5);
        Delta_T(j,i)=(Hax (j,i)*cos(I)*cos(A)+Hay(j,i)*cos(I)*sin(A)+Za(j,i)*sin(I));
    end
end

figure()
subplot 221; contourf(x,y,Hax,12); axis xy; xlabel({'X/m';'(a) Ha_x ��ֵ��ͼ'});ylabel('Y/m'); colorbar;
subplot 222; contourf(x,y,Hay,12); axis xy; xlabel({'X/m';'(b) Ha_y ��ֵ��ͼ'});ylabel('Y/m'); colorbar;
subplot 223; contourf(x,y,Za,12); axis xy; xlabel({'X/m';'(c) Z_a ��ֵ��ͼ'});ylabel('Y/m'); colorbar;
subplot 224; contourf(x,y,Delta_T,12); axis xy; xlabel({'X/m';'(d) \Delta_T ��ֵ��ͼ'});ylabel('Y/m'); colorbar;
