clc

% �ٶ���ʵ��ģ�Ͳ���
x0s=1000.0;
z0s=400.0;
ls=200;
bs=100;
alphas=pi/6;
sigmas=8.67;
XS=[x0s,z0s,ls,bs,alphas,sigmas];
% ģ�ͷ��ݳ�ʼ����
x0=1000.0;
z0=1000.0;
l=500;
b=500;
alpha=pi/2;
sigma=3;
X0=[x0,z0,l,b,alpha,sigma];
% ��������
N=101;      %�۲����
E=1e-7;     %������׼
NS=100;      %����������
% �۲������
XK=linspace(0,2000,N);
ZK=zeros(1,N);
	
% ʵ��ֵ
YS=zhengyan(XS,XK,ZK);
% ����
X1=fanyan(YS,X0,XK,ZK,E,NS)
