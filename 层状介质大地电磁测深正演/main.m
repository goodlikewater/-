clear all;
p = [100, 600];
h = 1800;
mu=(4e-7)*pi;
T=logspace(-3,4,40);
i=sqrt(-1);
k=zeros(size(p,2),size(T,2));
for N=1:size(p,2)
    k(N,:)=sqrt(-i*2*pi*mu./(T.*p(N)));
end
m=size (p,2);
y=-(i*mu*2*pi)./(T.*k(m,:));
%�ײ��迹
for n=m-1:-1:1
    A=-(i*2*pi*mu)./(T.*k(n,:));
    B=exp(-2*k(n,:)*h(n));
    y=A.*(A.*(1-B)+y.*(1+B))./(A.*(1+B)+y.*(1-B)); %���沨�迹
end
pc=(T./(mu*2*pi)).*(abs(y).^2); %�ӵ�����
ph=-atan(imag(y)./real(y)).*180/pi; %��λ

figure();
subplot 121;semilogx(T,pc);xlabel({'T/s';'(a) �ӵ���������'});ylabel('\rho_a/\Omega�� m');
subplot 122;semilogx(T,ph);xlabel({'T/s';'(b) ��λ����'});ylabel('Phase/(��)');
