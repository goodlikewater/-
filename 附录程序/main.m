clc

% 假定的实际模型参数
x0s=1000.0;
z0s=400.0;
ls=200;
bs=100;
alphas=pi/6;
sigmas=8.67;
XS=[x0s,z0s,ls,bs,alphas,sigmas];
% 模型反演初始参数
x0=1000.0;
z0=1000.0;
l=500;
b=500;
alpha=pi/2;
sigma=3;
X0=[x0,z0,l,b,alpha,sigma];
% 其他参数
N=101;      %观测点数
E=1e-7;     %收敛标准
NS=100;      %迭代最大次数
% 观测点坐标
XK=linspace(0,2000,N);
ZK=zeros(1,N);
	
% 实测值
YS=zhengyan(XS,XK,ZK);
% 反演
X1=fanyan(YS,X0,XK,ZK,E,NS)
