% 正演计算
function Y=zhengyan(X,XK,ZK)
 
% 参数
x0=X(1);
z0=X(2);
l=X(3);
b=X(4);
alpha=X(5);
sigma=X(6);
G=6.67e-5;
 
N=length(XK);
Y=1:N;
 
% 计算
z1=z0-l*sin(alpha);
z2=z0+l*sin(alpha);
for n=1:N
    xk=XK(n);
    zk=ZK(n);
    
    f1=pi-arctan(z1,xk-x0+b+l*cos(alpha));
    f2=pi-arctan(z2,xk-x0+b-l*cos(alpha));
    f3=pi-arctan(z1,xk-x0-b+l*cos(alpha));
    f4=pi-arctan(z2,xk-x0-b-l*cos(alpha));
    
    r1=sqrt((xk-x0+b+l*sin(alpha))^2+(z0-zk-l*sin(alpha))^2);
    r2=sqrt((xk-x0+b-l*sin(alpha))^2+(z0-zk+l*sin(alpha))^2);
    r3=sqrt((xk-x0-b+l*sin(alpha))^2+(z0-zk-l*sin(alpha))^2);
    r4=sqrt((xk-x0-b-l*sin(alpha))^2+(z0-zk+l*sin(alpha))^2);
    
    Y(n)=2*G*sigma*((z2*(f2-f4)-z1*(f1-f3))+xk*((sin(alpha))^2*log(r2*r3/(r1*r4))+cos(alpha)*sin(alpha)*(f1-f2-f3+f4))+2*b*((sin(alpha))^2*log(r4/r3)+cos(alpha)*sin(alpha)*(f3-f4)));
end
end
